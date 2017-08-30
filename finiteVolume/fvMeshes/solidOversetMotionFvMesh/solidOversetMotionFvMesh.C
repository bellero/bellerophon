/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "transformField.H"
#include "cellZoneMesh.H"
#include "boolList.H"
#include "syncTools.H"
#include "triSurfaceSearch.H"
#include "cellSet.H"
#include "faceSet.H"
#include "interpolation.H"
#include "PstreamBuffers.H"

#include "processorCyclicPolyPatch.H"

#include "solidOversetMotionFvMesh.H"
#include "bellerophonInterpolation.H"
#include "bellerophonInterface.H"
#include "bellerophonPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidOversetMotionFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, solidOversetMotionFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidOversetMotionFvMesh::solidOversetMotionFvMesh(const IOobject& io)
:
    bellerophonFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                this->time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::AUTO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    SBMFPtr_(solidBodyMotionFunction::New(dict_, io.time())),
    oldTransf_(SBMFPtr_().transformation()),
    points0_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    oversetPatchID_(this->boundary().findPatchID(dict_.lookup("oversetPatch"))),
    oversetZoneID_(-1),
    UName_(dict_.lookupOrDefault<word>("UName","U")),
    checkLiveCells_(dict_.lookupOrDefault<Switch>("checkLiveCells",false)),
    interpolateLiveCells_
    (
        dict_.lookupOrDefault<Switch>("interpolateLiveCells",false)
    ),
    interpolationSchemeName_
    (
        dict_.lookupOrDefault<word>("interpolationScheme","cellPoint")
    ),
    correctUf_
    (
        dict_.lookupOrDefault<Switch>("correctUf",true)
    ),
    UfName_(dict_.lookupOrDefault<word>("UfName","Uf")),
    interfaceIndex_(-1)
{
    if (points0_.size() != nPoints())
    {
        FatalIOErrorIn
        (
            "solidOversetMotionFvMesh::solidOversetMotionFvMesh(const IOobject&)",
            dict_
        )   << "Read " << points0_.size()
            << " undisplaced points from " << points0_.objectPath()
            << " but the current mesh has " << nPoints()
            << exit(FatalIOError);
    }

    if(!isA<bellerophonPolyPatch>(this->boundaryMesh()[oversetPatchID_]))
    {
        FatalErrorIn
        (
            "solidOversetMotionFvMesh::solidOversetMotionFvMesh(const IOobject&)"
        )
        << "Patch " << word(dict_.lookup("oversetPatch")) << " is not an "
        << "overset patch."
        << exit(FatalError);
    }

    const bellerophonPolyPatch& oversetPatch =
        refCast<const bellerophonPolyPatch>
        (
            this->boundaryMesh()[oversetPatchID_]
        );

    oversetZoneID_ = oversetPatch.oversetZoneID();

    // Recreate original mesh geometry and hole boundary
    {
        // Cache current point positions
        const pointField curPoints = this->points();

        // Move to original positions
        this->movePoints(points0_);

        // Generate hole boundary and save points
        holeBoundaryPtr_ = oversetPatch.holeBoundary();
        boundaryPoints0_ = holeBoundaryPtr_().points();

        if(debug)
        {
            holeBoundaryPtr_().write("initalPos.stl");
        }

        // Go back to current position
        this->movePoints(curPoints);

        // Update hole boundary
        holeBoundaryPtr_().movePoints
        (
            transformPoints
            (
                SBMFPtr_().transformation(),
                boundaryPoints0_
            )
        );
        // Possibly reset hole boundary instead of transformation
        //holeBoundaryPtr_.reset(oversetPatch.holeBoundary());
    }

    generateMotionMap();

    interfaceIndex_ = addHole();


    writeLiveCells(this->time().timeName());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidOversetMotionFvMesh::~solidOversetMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::labelList>
Foam::solidOversetMotionFvMesh::createHoleMap() const
{
    // Triangular surface search
    triSurfaceSearch searchEngine(holeBoundaryPtr_(), 1e-08, 10);

    // Is the cell inside the hole surface
    const boolList inHole( searchEngine.calcInside( this->cellCentres() ) );

    // Markup for hole cells
    autoPtr<labelList> mapPtr ( new labelList( this->nCells() ) );

    labelList& holeMap = mapPtr();

    // Mark hole cells
    forAll(holeMap, cellI)
    {
        if(inHole[cellI]) holeMap[cellI] = bellerophon::HOLE_CELL;
        else holeMap[cellI] = bellerophon::REGULAR_CELL;
    }

    // Remove moving cells
    forAll(movingCellIDs_, cellI)
    {
        holeMap[movingCellIDs_[cellI]] = bellerophon::REGULAR_CELL;
    }

    const labelList& own = this->owner();
    const labelList& nei = this->neighbour();

    forAll(own, faceI)
    {
        register const label o = own[faceI];
        register const label n = nei[faceI];
        if
        (
            holeMap[o] == bellerophon::HOLE_CELL
            &&
            holeMap[n] == bellerophon::REGULAR_CELL
        )
        {
            holeMap[o] = bellerophon::ACCEPTOR_CELL;
        }
        else if
        (
            holeMap[o] == bellerophon::REGULAR_CELL
            &&
            holeMap[n] == bellerophon::HOLE_CELL
        )
        {
            holeMap[n] = bellerophon::ACCEPTOR_CELL;
        }
    }

    // Account for processor patches
    // i.e. mark hole cells as acceptor cells, if cell behind processor
    // interface is regular cell
    if(Pstream::parRun())
    {
        const label nIntFaces = this->nInternalFaces();

        // Copy state to separate list to swap
        labelList neighbourCell
        (
            this->nFaces()-nIntFaces,
            bellerophon::REGULAR_CELL
        );

        const polyBoundaryMesh& patches = this->boundaryMesh();

        forAll(patches,patchI)
        {
            const polyPatch& patch = patches[patchI];
            if
            (
                isA<processorPolyPatch>(patch)
                &&
                !isA<processorCyclicPolyPatch>(patch)
            )
            {
                const labelList fcs = patch.faceCells();

                const label start = patch.start() - nIntFaces;
                forAll(fcs, faceI)
                {
                    neighbourCell[start+faceI] = holeMap[fcs[faceI]];
                }
            }
        }

        // Swap state
        syncTools::swapBoundaryFaceList(*this, neighbourCell);

        // If only side marks boundary face it is part of the interface,
        // otherwise not
        forAll(patches,patchI)
        {
            const polyPatch& patch = patches[patchI];
            if
            (
                isA<processorPolyPatch>(patch)
                &&
                !isA<processorCyclicPolyPatch>(patch)
            )
            {
                const labelList fcs = patch.faceCells();

                const label start = patch.start() - nIntFaces;
                forAll(fcs, faceI)
                {
                    if
                    (
                        neighbourCell[start+faceI] == bellerophon::REGULAR_CELL
                        &&
                        holeMap[fcs[faceI]] == bellerophon::HOLE_CELL
                    )
                    holeMap[fcs[faceI]] = bellerophon::ACCEPTOR_CELL;
                }
            }
        }
    }

    return mapPtr;
}

void Foam::solidOversetMotionFvMesh::generateMotionMap()
{
    // Moving cells mask
    boolList moveCells(this->nCells(), false);

    label nMovingCells = 0;

    // Cells in moving zone
    const labelList& cellIDs = this->cellZones()[oversetZoneID_];

    forAll(cellIDs, i)
    {
        moveCells[cellIDs[i]] = true;
        nMovingCells++;
    }

    movingCellIDs_.setSize(nMovingCells);
    nMovingCells = 0;

    // Moving points mask
    boolList movePts(this->nPoints(), false);

    forAll(moveCells, cellI)
    {
        if(moveCells[cellI])
        {
            const cell& c = this->cells()[cellI];
            forAll(c, j)
            {
                const face& f = this->faces()[c[j]];
                forAll(f, k)
                {
                    label pointI = f[k];
                    movePts[pointI] = true;
                }
            }
            movingCellIDs_[nMovingCells++] = cellI;
        }
    }

    moveCells.setSize(0);

    syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

    // Count moving points
    label nMovingPoints = 0;
    forAll(movePts, i)
    {
        if (movePts[i])
        {
            nMovingPoints++;
        }
    }

    movingPointIDs_.setSize(nMovingPoints);
    nMovingPoints=0;

    // File moving point list
    forAll(movePts, i)
    {
        if (movePts[i])
        {
            movingPointIDs_[nMovingPoints++]=i;
        }
    }
}


Foam::label Foam::solidOversetMotionFvMesh::addHole()
{
    autoPtr<labelList> mapPtr = createHoleMap();
    const labelList& holeMap = mapPtr();

    bellerophon::Interpolation().markHole();

    bellerophon::Interpolation().clearTopology();

    return bellerophon::Interpolation().addHole
    (
        holeMap,
        oversetZoneID_,
        oversetPatchID_
    );
}


void Foam::solidOversetMotionFvMesh::interpolateNewLiveCells
(
    const boolList& map,
    const pointField& transformedPos
)
{
    // Number of cells to interpolate
    label centreI = 0;
    forAll(map, cellI) if(map[cellI]) centreI++;

    // Number of this processor
    const label myProcNo = Pstream::myProcNo();

    // Number of this processor
    const label nProcs = Pstream::nProcs();

    autoPtr< List<searchItem> > itemsPtr(new List<searchItem>(centreI));
    {
        List<searchItem>& items = itemsPtr();

        centreI = 0;
        forAll(map, cellI)
        {
            if(map[cellI])
            {
                searchItem& item = items[centreI];
                item.procID() = myProcNo;
                item.cellLabel() = cellI;
                item.position() = transformedPos[centreI];
                item.seed() = -1;
                item.zoneID() = oversetZoneID_;
                centreI++;

            }
        }
    }

    gradientSearch gs(*this);

    gs.search(itemsPtr);

    labelList nItemsPerProc(nProcs,0);

    List<searchItem>& items = itemsPtr();

    forAll(items, itemI)
    {
        nItemsPerProc[items[itemI].procID()]++;
    }

    labelListList IDperProc(nProcs);
    forAll(IDperProc, procI)
    {
        IDperProc[procI].setSize(nItemsPerProc[procI]);
        nItemsPerProc[procI] = 0;
    }

    forAll(items, itemI)
    {
        const label procI = items[itemI].procID();
        IDperProc[procI][nItemsPerProc[procI]++] = itemI;
    }

    // Interpolate all types
    interpolateFields<scalar>(items, IDperProc);
    interpolateFields<vector>(items, IDperProc);
    interpolateFields<tensor>(items, IDperProc);
    interpolateFields<symmTensor>(items, IDperProc);
    interpolateFields<sphericalTensor>(items, IDperProc);

    itemsPtr.clear();
}

template<class Type>
void Foam::solidOversetMotionFvMesh::interpolateFields
(
    const List<searchItem>& items,
    const labelListList& IDtoProc
)
{
    const label myProcNo = Pstream::myProcNo();
    const label nProcs = Pstream::nProcs();

    // For all fields of that type
    const wordList fieldNames =
        this->names< GeometricField<Type, fvPatchField, volMesh> >();


    forAll(fieldNames, nameI)
    {
        if(debug)
        {
            Info<<"Interpolating field "<<fieldNames[nameI]<<endl;
        }

        const GeometricField<Type, fvPatchField, volMesh>& field =
            this->lookupObject<GeometricField<Type, fvPatchField, volMesh> >
            (
                fieldNames[nameI]
            );

        List<Field<Type> > interpolated(nProcs);

        autoPtr<interpolation<Type> > interpolatorPtr =
            interpolation<Type>::New(interpolationSchemeName_, field);

        interpolation<Type>& interpolator = interpolatorPtr();

        labelListList itemToCell(nProcs);

        forAll(IDtoProc, procI)
        {
            itemToCell[procI].setSize(IDtoProc[procI].size());
            interpolated[procI].setSize(IDtoProc[procI].size());
        }

        labelList nItemsPerProc(nProcs, 0);

        forAll(items, itemI)
        {
            const searchItem& item = items[itemI];
            const label procI = item.procID();

            itemToCell[procI][nItemsPerProc[procI]] = item.cellLabel();

            interpolated[procI][nItemsPerProc[procI]++] =
                interpolator.interpolate(item.position(), item.seed());
        }

        // Send interpolated and asign to internal field
        PstreamBuffers pBufs(Pstream::nonBlocking);

        forAll(interpolated, procI)
        {
            if(procI != myProcNo)
            {
                UOPstream toDomain(procI, pBufs);
                toDomain << itemToCell[procI] << interpolated[procI];
            }
        }

        pBufs.finishedSends();

        forAll(interpolated, procI)
        {
            if(procI != myProcNo)
            {
                UIPstream str(procI, pBufs);
                str >> itemToCell[procI] >> interpolated[procI];
            }
        }

        Field<Type>& iField =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
                (field).primitiveFieldRef();

        forAll(itemToCell, procI)
        {
            const labelList& toCell = itemToCell[procI];
            const Field<Type>& values = interpolated[procI];
            forAll(toCell, cellI)
            {
                iField[toCell[cellI]] = values[cellI];
            }
        }
    }
}

bool Foam::solidOversetMotionFvMesh::update()
{
    // Transformation
    const septernion transf = SBMFPtr_().transformation();

    // New mesh vertex positions
    pointField newPoints(points0_);

    // autoPtr to old hole cells
    autoPtr<boolList> oldHoleCellsPtr;

    // Old acceptor cells
    const labelList oldAcceptors(bellerophon::Interpolation().acceptorCells());

    if(checkLiveCells_)
    {
        // Old hole cells
        oldHoleCellsPtr.reset(new boolList(this->nCells(), false));

        boolList& map = oldHoleCellsPtr();

        // Reference to hole cells, changes when updating hole
        const labelList& holeCells = bellerophon::Interpolation().holeCells();

        forAll(holeCells, cellI) map[holeCells[cellI]] = true;
    }

    // Transform
    UIndirectList<point>(newPoints, movingPointIDs_) =
        transformPoints(transf, pointField(newPoints, movingPointIDs_));

    // Move mesh
    fvMesh::movePoints(newPoints);

    // Update position of hole boundary surface
    holeBoundaryPtr_().movePoints(transformPoints(transf, boundaryPoints0_));

    interfaceIndex_ = addHole();

    // List for internal faces next to cells that have been hole cells and are
    // now next to live cells
//     labelList interpolateFaces(0);
    if(checkLiveCells_)
    {
        boolList& map = oldHoleCellsPtr();

        // Reference to hole cells, changes when updating hole
        const labelList& holeCells = bellerophon::Interpolation().holeCells();

        const labelList& acceptorCells =
        bellerophon::Interpolation().acceptorCells();

        forAll(holeCells, cellI) map[holeCells[cellI]] = false;

        forAll(acceptorCells, acceptorI) map[acceptorCells[acceptorI]] = false;

        // TODO consider processor boundary faces

        label globalNewLiveCells = 0;

        forAll(map, cellI) if(map[cellI]) globalNewLiveCells++;

        pointField transformedPos(globalNewLiveCells);

        const label localNewLiveCells = globalNewLiveCells;

        globalNewLiveCells = 0;

        const pointField& ccs = this->cellCentres();

        forAll(map, cellI)
        {
            if(map[cellI]) transformedPos[globalNewLiveCells++] = ccs[cellI];
        }

        // Transform positions of new live cells to positions of moved overset
        // grid cells
        transformedPos = transformPoints(inv(oldTransf_),transformedPos);
        transformedPos = transformPoints(transf,transformedPos);

        Pstream::gather<label>(globalNewLiveCells, sumOp<label>());
        Pstream::scatter<label>(globalNewLiveCells);

        Info<<"    There are "<<globalNewLiveCells<<" new live cells out of "
            <<holeCells.size()<<", that used to be hole cells."<<endl;

        if(globalNewLiveCells > 0 && interpolateLiveCells_)
        {
            Info<<"    Interpolating from old overset cells."<<endl;

            interpolateNewLiveCells(map, transformedPos);
        }
    }

    if( (debug && this->time().outputTime()) || debug>1)
    {
        holeBoundaryPtr_().write
        (
            "hole"+Foam::name(Pstream::myProcNo())+"_"
           +this->time().timeName()+".stl"
        );
    }

    static bool hasWarned = false;

    if (foundObject<volVectorField>(UName_))
    {
        const volVectorField& U = lookupObject<volVectorField>(UName_);

        // Correct surface velocity field on faces between live cells and
        // acceptor cells where
        if
        (
            correctUf_
        )
        {
            if(!foundObject<surfaceVectorField>(UfName_))
            {
                WarningIn
                (
                    "bool Foam::solidOversetMotionFvMesh::update()"
                )
                << "Cannot find surfaceVectorField " << UfName_
                << " to correct face velocities."  << endl;
            }
            else
            {
                // Interpolate velocities to acceptor cells
                bellerophon::Interpolation().interpolate<vector>
                (
                    UName_,
                    Pstream::nonBlocking,
                    false
                );

                // Interpolate velocities to acceptor cells
                bellerophon::Interpolation().interpolate<vector>
                (
                    UName_+"_0",
                    Pstream::nonBlocking,
                    false
                );

                // Get new acceptor cells
                boolList map(this->nCells(), false);

                // Interpolate p to new acceptor cells to get gradient for
                // momentum equations right
                {
                    volScalarField& p = const_cast<volScalarField&>
                        (
                            lookupObject<volScalarField>("p")
                        );

                    volScalarField interpolatedP("interpolatedP", p);
                    interpolatedP.checkIn();

                    // Interpolate pressure values to acceptor cells
                    bellerophon::Interpolation().interpolate<scalar>
                    (
                        "interpolatedP",
                        Pstream::nonBlocking,
                        false
                    );

                    const labelList& newAcceptors =
                        bellerophon::Interpolation().acceptorCells();

                    forAll(newAcceptors, cellI)
                    {
                        map[newAcceptors[cellI]] = true;
                    }

                    forAll(oldAcceptors, cellI)
                    {
                        map[oldAcceptors[cellI]] = false;
                    }

                    forAll(map, cellI)
                    {
                        if(map[cellI])
                        {
                            p.primitiveFieldRef()[cellI] =
                                interpolatedP.internalField()[cellI];
                        }
                    }
                }

                // Faces between live cells and acceptors
                const labelList& interfaceFaces =
                    bellerophon::Interpolation().interfaceFaces
                    (
                        interfaceIndex_
                    );

                // Flip map for faces between live cells and acceptors
                const boolList& interfaceFlipMap =
                    bellerophon::Interpolation().interfaceFlipMap
                    (
                        interfaceIndex_
                    );

                Info<<"    Correcting "<<UfName_<<" on "<<interfaceFaces.size()
                    <<" faces."<<endl;


                // Access to surface velocities
                surfaceVectorField& Uf =
                    const_cast<surfaceVectorField&>
                    (
                        lookupObject<surfaceVectorField>(UfName_)
                    );

                surfaceVectorField& UfOld = Uf.oldTime();

                const surfaceVectorField UoldF = fvc::interpolate(U.oldTime());
                const vectorField& intUoldF = UoldF.internalField();

                // Surface velocities on internal faces
                vectorField& intUf = Uf.primitiveFieldRef();
                vectorField& intUfOld = UfOld.primitiveFieldRef();

                const labelList& own = this->owner();
                const labelList& nei = this->neighbour();

                // Face fluxes from face velocities
                tmp<surfaceScalarField> tPhi =
                    this->Sf() & Uf;

                // Divergence from calculated face fluxes
                tmp<volScalarField> tDivPhi = fvc::div(tPhi());

                // Cell values of divergence of face fluxes
                const scalarField& intDivPhi = tDivPhi().internalField();

                // Sizes of faces
                const vectorField& intSf = this->Sf().internalField();

                // Sizes of faces
                const scalarField& intMagSf =
                    this->magSf().internalField();

                const scalarField& cellVols = this->V();

                // Sum up interface face sums for normalization
                scalarField faceSums(own.size(), 0.0);

                boolList adjustMap(own.size(), false);
                boolList flipMap(own.size(), false);

                forAll(interfaceFaces, faceI)
                {
                    const label f = interfaceFaces[faceI];
                    adjustMap[f] = true;
                    flipMap[f] = interfaceFlipMap[faceI];
                }

                //TODO
#warning does adjustMap has to be set to false here???

                forAll(own, faceI)
                {
                    if(map[own[faceI]] == map[nei[faceI]])
                    {
                        flipMap[faceI] = false;
                    }
                }

                forAll(adjustMap, faceI)
                {
                    if(adjustMap[faceI])
                    {
                        const bool   flip  = flipMap[faceI];
                        const label  cellI = flip ? nei[faceI] : own[faceI];
                        faceSums[cellI] += intMagSf[faceI];
                    }
                }

                forAll(adjustMap, faceI)
                {
                    if(adjustMap[faceI])
                    {
                        const vector n      = intSf[faceI] / intMagSf[faceI];
                        const bool   flip   = flipMap[faceI];
                        const label  cellI  = flip ? nei[faceI] : own[faceI];
                        const scalar vol    = cellVols[cellI];
                        const scalar area   =
                            flip ? faceSums[cellI] : -faceSums[cellI];
                        const vector factor = vol * n / area;

                        intUf[faceI] += intDivPhi[cellI] * factor;
                        intUfOld[faceI] = intUoldF[faceI];
                    }
                }
//                 forAll(interpolateFaces, faceI)
//                 {
//                     const label f = interpolateFaces[faceI];
//                     internalUf[f] = internalNewSf[f];
//                 }

            }
        }

        // Now correct boundary conditions of velocity field
        const_cast<volVectorField&>(U).correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningIn("solidOversetMotionFvMesh::update()")
            << "Could not find volVectorField " << UName_
            << ". Therefore not going to updating boundary conditions."
            << endl;
    }

    oldTransf_ = transf;
    return true;
}


bool Foam::solidOversetMotionFvMesh::writeLiveCells(const word& timeName) const
{
    const labelList& holeCells = bellerophon::Interpolation().holeCells();

    const labelList& acceptorCells =
    bellerophon::Interpolation().acceptorCells();

    // Count live cells
    label nLiveCells = this->nCells();

    // Markup for live cells
    boolList liveCellMap(nLiveCells, true);

    // Remove hole cells from live cells
    forAll(holeCells,cellI)
    {
        if(liveCellMap[holeCells[cellI]] == true)
        {
            liveCellMap[holeCells[cellI]] = false;
            nLiveCells--;
        }
    }

    // Cell set for output
    cellSet liveCellSet
    (
        *this,
        "liveCells",
        nLiveCells
    );

    cellSet acceptorCellSet
    (
        *this,
        "acceptorCells",
        acceptorCells.size()
    );

    // Remove acceptor cells from live cells
    forAll(acceptorCells,cellI)
    {
        acceptorCellSet.insert(acceptorCells[cellI]);
        if(liveCellMap[acceptorCells[cellI]] == true)
        {
            liveCellMap[acceptorCells[cellI]] = false;
            nLiveCells--;
        }
    }

    // Add cells to cell set
    forAll(liveCellMap, cellI) if(liveCellMap[cellI]) liveCellSet.insert(cellI);

    // Set instance to current time
    liveCellSet.instance() = timeName;
    acceptorCellSet.instance() = timeName;

    return liveCellSet.write() && acceptorCellSet.write();
}


bool Foam::solidOversetMotionFvMesh::writeObject
(
    Foam::IOstream::streamFormat fmt,
 Foam::IOstream::versionNumber ver,
 Foam::IOstream::compressionType cmp
) const
{
    return writeLiveCells(this->time().timeName()) && this->dynamicFvMesh::writeObject(fmt, ver, cmp);
}

// ************************************************************************* //
