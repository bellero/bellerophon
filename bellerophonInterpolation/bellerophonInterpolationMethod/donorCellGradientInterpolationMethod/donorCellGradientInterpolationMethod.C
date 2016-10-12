/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "donorCellGradientInterpolationMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "searchItem.H"
#include "PstreamBuffers.H"
#include "surfaceMesh.H"
#include "processorCyclicFvPatch.H"
#include "fvsPatchFields.H"

namespace Foam
{
    defineTypeNameAndDebug(donorCellGradientInterpolationMethod, 0);
    addToRunTimeSelectionTable(bellerophonInterpolationMethod,donorCellGradientInterpolationMethod, bellerophonInterpolationMethod);
}

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

Foam::donorCellGradientInterpolationMethod::donorCellGradientInterpolationMethod
(
    const Foam::fvMesh& mesh,
    const Foam::labelList& primaryDonorCells,
    const Foam::labelList& acceptorCells,
    const Foam::vectorField& deltas,
    const Foam::autoPtr<Foam::List< Foam::searchItem> >& donorItemsPtr,
    const dictionary& dict
)
:
bellerophonInterpolationMethod
(
    mesh,
    primaryDonorCells,
    acceptorCells,
    deltas,
    donorItemsPtr,
    dict
)
{}

// * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * //

Foam::donorCellGradientInterpolationMethod::
~donorCellGradientInterpolationMethod()
{}

// * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

void Foam::donorCellGradientInterpolationMethod::update()
{
    primaryDonorWeightsPtr_.reset(new scalarField(0));

    // Search extended donor cells if using implicit interpolation
    // Count interpolation cells on this proc and on neighbour procs
    label nOwnProcInterpolationCells = 0;

    // List of faces for each cell
    const cellList& cells = mesh_.cells();

    // Number of internal faces
    const label nIntFaces = mesh_.nInternalFaces();

    // Owner of faces
    const labelList& own = mesh_.owner();

    // Neighbour of faces
    const labelList& nei = mesh_.neighbour();

    const bool parRun = Pstream::parRun();

    const scalar nProcs = Pstream::nProcs();

    // Map boundary face to processor if is on processorPatch
    labelList faceToProcMap(0);
    labelList procToPatchMap(nProcs,-1);

    if(parRun)
    {
        faceToProcMap.setSize(mesh_.nFaces() - nIntFaces, -1);

        // Count processor patches with faces
        forAll(mesh_.boundary(), patchI)
        {
            if
            (
                mesh_.boundary()[patchI].size() > 0
                &&
                isA<processorFvPatch>(mesh_.boundary()[patchI])
                &&
                !isA<processorCyclicFvPatch>(mesh_.boundary()[patchI])
            )
            {
                label faceI = mesh_.boundary()[patchI].start() - nIntFaces;
                const label lastFace =
                    faceI + mesh_.boundary()[patchI].size();
                const label procI = refCast<const processorFvPatch>
                    (mesh_.boundary()[patchI]).neighbProcNo();

                while(faceI < lastFace)
                {
                    faceToProcMap[faceI++]=procI;
                }

                procToPatchMap[procI] = patchI;
            }
        }
    }

    // Count faces per proc
    labelList nNeighbourProcInterpolationCells(nProcs,0);
    forAll(primaryDonorCells_,donorI)
    {
        const label curDonor = primaryDonorCells_[donorI];
        const labelList& cellFaces = cells[curDonor];
        forAll(cellFaces, faceI)
        {
            const label curFace = cellFaces[faceI];
            if(curFace < nIntFaces)
            {
                nOwnProcInterpolationCells++;
            }
            else if
            (
                parRun
                &&
                faceToProcMap[curFace-nIntFaces] > -1
            )
            {
                nNeighbourProcInterpolationCells
                    [faceToProcMap[curFace-nIntFaces]]++;
            }
        }
    }

    ownInterpolationItemsPtr_.reset
    (
        new List<interpolationItem>(nOwnProcInterpolationCells)
    );

    neighbourInterpolationItemsPtr_.reset
    (
        new List< List<interpolationItem> >(nProcs)
    );

    neighbourValueToFieldMapPtr_.reset
    (
        new labelListList(nProcs)
    );

    scalarField& primaryDonorWeights =
        primaryDonorWeightsPtr_();

    List<interpolationItem>& ownInterpolationItems =
        ownInterpolationItemsPtr_();

    List< List<interpolationItem> >& neighbourInterpolationItems =
        neighbourInterpolationItemsPtr_();

    // Maps interpolation contribution of interpolation items on
    // neighbour proc to interpolated field
    labelListList& neighbourValueToFieldMap =
        neighbourValueToFieldMapPtr_();

    const label nDonors = primaryDonorCells_.size();
    primaryDonorWeights.setSize(nDonors);

    nOwnProcInterpolationCells = 0;

    forAll(nNeighbourProcInterpolationCells, procI)
    {
        neighbourInterpolationItems[procI].setSize
        (
            nNeighbourProcInterpolationCells[procI]
        );
        neighbourValueToFieldMap[procI].setSize
        (
            nNeighbourProcInterpolationCells[procI]
        );
        nNeighbourProcInterpolationCells[procI] = 0;
    }

    const surfaceScalarField& faceWeights =
        mesh_.weights();

    const scalarField& internalFaceWeights =
        faceWeights.internalField();

    const surfaceVectorField& Sf = mesh_.Sf();

    const vectorField& internalSf = Sf.internalField();

    const scalarField& V = mesh_.V().field();

    forAll(primaryDonorCells_,donorI)
    {
        const label curDonor = primaryDonorCells_[donorI];

        // Interpolation item for primary donor cell
        scalar& primaryDonorWeight =
            primaryDonorWeights[donorI];

        primaryDonorWeight = 1.0;

        // Const factor for weight calculation
        const vector weightVector = deltas_[donorI] / V[curDonor];

        const labelList& cellFaces = cells[curDonor];
        forAll(cellFaces, faceI)
        {
            const label curFace = cellFaces[faceI];

            if(curFace < nIntFaces)
            {
                // Interpolation item for neighbour cell
                interpolationItem& neighbourItem =
                    ownInterpolationItems[nOwnProcInterpolationCells++];

                // Weight factor for the surface value
                // i.e. Delta * 1/Vp * Sf
                const scalar weightFactor =
                    weightVector & internalSf[curFace];

                // Interpolation weight for the face value
                // i.e. lambda
                const scalar faceWeight = internalFaceWeights[curFace];

                if(curDonor == own[curFace])
                {
                    // Primary donor cell is owner of face
                    primaryDonorWeight += faceWeight * weightFactor;
                    neighbourItem.set
                    (
                        donorI,
                        nei[curFace],
                        (1.0 - faceWeight) * weightFactor
                    );
                }
                else
                {
                    // Primary donor cell is neighbour of face
                    primaryDonorWeight +=
                        (faceWeight - 1.0) * weightFactor;
                    neighbourItem.set
                    (
                        donorI,
                        own[curFace],
                        - faceWeight * weightFactor
                    );
                }
            }
            else if
            (
                parRun
                &&
                faceToProcMap[curFace-nIntFaces] > -1
            )
            {
                const label procI = faceToProcMap[curFace-nIntFaces];

                // Interpolation item for neighbour cell
                interpolationItem& neighbourItem =
                    neighbourInterpolationItems[procI]
                    [nNeighbourProcInterpolationCells[procI]];

                // Label of the patch
                const label patchI =
                    mesh_.boundaryMesh().patchID()[curFace-nIntFaces];

                // Number of the face on the patch
                const label patchFace =
                    curFace - mesh_.boundaryMesh()[patchI].start();

                // Weight factor for the surface value
                // i.e. Delta * 1/Vp * Sf
                const scalar weightFactor = weightVector &
                    Sf.boundaryField()[patchI][patchFace];

                // Interpolation weight for the face value
                // i.e. lambda
                const scalar faceWeight =
                    faceWeights.boundaryField()[patchI][patchFace];

                primaryDonorWeight += faceWeight * weightFactor;
                neighbourItem.set
                (
                    donorI,
                    patchFace, // setting face on the patch as donor cell
                                // will be converted to cell behind this
                                // face when received by neighbour proc
                    (1.0 - faceWeight) * weightFactor
                );
                neighbourValueToFieldMap[procI]
                    [nNeighbourProcInterpolationCells[procI]] =
                    donorI;
                nNeighbourProcInterpolationCells[procI]++;
            }
            else
            {
                // Treat every patch a zeroGradient for now, fell free to add
                // others to RHS

                // Label of the patch
                const label patchI =
                    mesh_.boundaryMesh().patchID()[curFace-nIntFaces];

                if(!isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchI]))
                {
                    // Number of the face on the patch
                    const label patchFace =
                        curFace - mesh_.boundaryMesh()[patchI].start();

                    // Weight factor for the surface value
                    // i.e. Delta * 1/Vp * Sf
                    const scalar weightFactor = weightVector &
                        Sf.boundaryField()[patchI][patchFace];

                    // Face weight is 1 for zero gradient
                    primaryDonorWeight += weightFactor;
                }
            }
        }
    }

    distribute(procToPatchMap);
}

