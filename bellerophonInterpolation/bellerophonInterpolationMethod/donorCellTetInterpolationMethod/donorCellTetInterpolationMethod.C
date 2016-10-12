/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "donorCellTetInterpolationMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "processorCyclicFvPatch.H"

namespace Foam
{
    defineTypeNameAndDebug(donorCellTetInterpolationMethod, 0);
    addToRunTimeSelectionTable(bellerophonInterpolationMethod,donorCellTetInterpolationMethod, bellerophonInterpolationMethod);
}

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

Foam::donorCellTetInterpolationMethod::donorCellTetInterpolationMethod
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

Foam::donorCellTetInterpolationMethod::
~donorCellTetInterpolationMethod()
{}

// * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

void Foam::donorCellTetInterpolationMethod::update()
{
    // First check in which cell centre - face edge ends - face centre tet
    // the cell point is
    // If tet is found, check side of cell centre - face centre - edge midpoint
    // triangle

    // Cell centres
    const vectorField& cc = mesh_.C().internalField();

    // Face centres
    const vectorField& fc = mesh_.faceCentres();

    // Cells of the mesh
    const cellList& cells = mesh_.cells();

    // Faces of the mesh
    const faceList& faces = mesh_.faces();

    // Edges of the faces
    const edgeList& edges = mesh_.edges();
    // Points of the mesh
    const pointField& points = mesh_.points();

    // Edges per face
    const labelListList& faceEdges = mesh_.faceEdges();

    // Adjoint cells per point
    const labelListList& pointFaces = mesh_.pointFaces();

    const label nIntFaces = mesh_.nInternalFaces();

    const label nDonors = primaryDonorCells_.size();

    primaryDonorWeightsPtr_.reset(new scalarField(nDonors));

    scalarField& primaryWeights = primaryDonorWeightsPtr_();

    // Mesh point that marks corner of the tet containing the face
    labelList tetPoints(nDonors,-1);

    const bool parRun = Pstream::parRun();

    const scalar nProcs = Pstream::nProcs();

    const labelList& own = mesh_.owner();
    const labelList& nei = mesh_.neighbour();

    // Map boundary face to processor if is on processorPatch
    // TODO: store this a member data and update on demand?
    labelList faceToProcMap(mesh_.nFaces() - nIntFaces, -1);

    // Map processor number to boundary patch id
    labelList procToPatchMap(nProcs,-1);

    // Number of neighbour items per proc
    labelList nNeighbourItemsPerProc(nProcs,0);

    // Number of items per donor
    labelList nItemsPerDonor(nDonors,0);

    // Number of items per donor
    labelList nNeighbourItemsPerDonor(nDonors,0);

    if(parRun)
    {
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

    // Number of interpolation items
    //  (make an educated guess by summing um neighbour cells)
    label nOwnInterpolationItems = 0;

    // Map boundary face to patch
    const labelList& faceToPatch = mesh_.boundaryMesh().patchID();

    forAll(primaryDonorCells_, donorI)
    {
        const label curDonor = primaryDonorCells_[donorI];
        const cell& curCell = cells[curDonor];

        // Centre of the current donor cell
        const vector& curCentre = cc[curDonor];

        // Vector from the current donor cell to the acceptor cell
        const vector& delta = deltas_[donorI];

        primaryWeights[donorI] = 1.0/(SMALL+mag(delta));

        bool notFound = true;

        for(label faceI = 0; faceI<curCell.size() && notFound; faceI++)
        {
            const label curFaceI = curCell[faceI];

            const vector centreFaceVec = fc[curFaceI] - curCentre;

            const face curFace = faces[curFaceI];
            const labelList& curFaceEdges = faceEdges[curFaceI];
            for
            (
                label edgeI = 0;
                edgeI<curFaceEdges.size() && notFound;
                edgeI++
            )
            {
                const edge& curEdge = edges[curFaceEdges[edgeI]];

                const vector& start = points[curEdge.start()];
                const vector& end = points[curEdge.end()];

                const vector centreStartVec = start - curCentre;
                const vector centreEndVec = end - curCentre;

                const vector normalCSF = centreStartVec ^ centreFaceVec;
                const vector normalCFE = centreFaceVec ^ centreEndVec;
                const vector normalCES = centreEndVec ^ centreStartVec;

                bool flip = (centreStartVec & normalCFE) > 0.0;
                if
                (
                    ((normalCSF & delta) >= 0.0) == flip
                    &&
                    ((normalCFE & delta) >= 0.0) == flip
                    &&
                    ((normalCES & delta) >= 0.0) == flip
                )
                {
                    // Inside face tri tet, check side
                    notFound = false;

                    const vector midPoint = 0.5*(end - start);
                    const vector centreMidVec = midPoint - curCentre;
                    const vector normalCFM = centreFaceVec ^ centreMidVec;

                    label pointI;
                    if(((normalCFM & delta) >= 0.0) == flip)
                    {
                        // In start side tet
                        pointI = curEdge.start();
                    }
                    else
                    {
                        // In end side tet
                        pointI = curEdge.end();
                    }
                    tetPoints[donorI] = pointI;

                    const labelList& curPointFaces = pointFaces[pointI];

                    forAll(curPointFaces,pointFaceI)
                    {
                        const label curPointFace = curPointFaces[pointFaceI];
                        if
                        (
                            curPointFace < nIntFaces
                        )
                        {
                            if
                            (
                                own[curPointFace] == curDonor
                                ||
                                nei[curPointFace] == curDonor
                            )
                            {
                                nOwnInterpolationItems++;
                                nItemsPerDonor[donorI]++;
                            }
                        }
                        else
                        {
                            const label bFaceI = curPointFace-nIntFaces;

                            if
                            (
                                faceToProcMap[bFaceI] != -1
                            )
                            {
                                const fvPatch& p =
                                    mesh_.boundary()[faceToPatch[bFaceI]];

                                if
                                (
                                    p.faceCells()[curPointFace - p.start()]
                                    ==
                                    curDonor
                                )
                                {
                                    nNeighbourItemsPerProc[faceToProcMap[bFaceI]]++;
                                    nNeighbourItemsPerDonor[donorI]++;
                                }
                            }
                        }
                    }
                }
            }
            if(notFound)
            {
                tetPoints[donorI] = -1;
            }
        }
    }

    ownInterpolationItemsPtr_.reset
    (
        new List<interpolationItem>(nOwnInterpolationItems)
    );
    List<interpolationItem>& ownItems = ownInterpolationItemsPtr_();

    neighbourInterpolationItemsPtr_.reset
    (
        new List< List<interpolationItem> >(Pstream::nProcs())
    );
    List< List<interpolationItem> >& neighbourItems =
        neighbourInterpolationItemsPtr_();

    neighbourValueToFieldMapPtr_.reset
    (
        new labelListList(nProcs)
    );

    // Maps interpolation contribution of interpolation items on
    // neighbour proc to interpolated field
    labelListList& neighbourValueToFieldMap =
        neighbourValueToFieldMapPtr_();    nOwnInterpolationItems = 0;

    forAll(nNeighbourItemsPerProc, procI)
    {
        label& nNeighbourItems = nNeighbourItemsPerProc[procI];
        neighbourItems[procI].setSize(nNeighbourItems);
        neighbourValueToFieldMap[procI].setSize(nNeighbourItems);
        nNeighbourItems = 0;
    }

    forAll(primaryDonorCells_, donorI)
    {
        const label curDonor = primaryDonorCells_[donorI];

        if(tetPoints[donorI] != -1)
        {
            const labelList& curPointFaces = pointFaces[tetPoints[donorI]];

            scalar weightsum = primaryWeights[donorI];

            label nItems = nItemsPerDonor[donorI];

            scalarList weights(nItems);
            labelList donors(nItems);

            label nNeighbourItems = nNeighbourItemsPerDonor[donorI];

            scalarList neighbourWeights(nNeighbourItems);
            labelList neighbourFace(nNeighbourItems);
            labelList neighbourProc(nNeighbourItems);

            const vector& acceptorPos = donorItemsPtr_()[donorI].position();

            forAll(curPointFaces,pointFaceI)
            {
                const label curPointFace = curPointFaces[pointFaceI];
                if
                (
                    curPointFace < nIntFaces
                )
                {
                    if
                    (
                        own[curPointFace] == curDonor
                    )
                    {
                        nItems--;
                        donors[nItems] = nei[curPointFace];
                        weights[nItems] =
                            1.0/(SMALL + mag( cc[nei[curPointFace]]-acceptorPos) );
                        weightsum += weights[nItems];
                    }
                    else if
                    (
                        nei[curPointFace] == curDonor
                    )
                    {
                        nItems--;
                        donors[nItems] = own[curPointFace];
                        register const scalar w = 1.0 /
                            (SMALL + mag( cc[own[curPointFace]]-acceptorPos) );
                        weights[nItems] = w;
                        weightsum += w;
                    }
                }
                else
                {
                    const label bFaceI = curPointFace-nIntFaces;

                    if
                    (
                        faceToProcMap[bFaceI] != -1
                    )
                    {
                        const fvPatch& p =
                            mesh_.boundary()[faceToPatch[bFaceI]];

                        if
                        (
                            p.faceCells()[curPointFace - p.start()]
                            ==
                            curDonor
                        )
                        {
                            nNeighbourItems--;

                            const label procID = faceToProcMap[bFaceI];

                            const processorPolyPatch& neiPatch =
                                refCast<const processorFvPatch>
                                    (
                                        mesh_.boundary()[procToPatchMap[procID]]
                                    ).procPolyPatch();

                            const label neighbourFaceI =
                                curPointFace-neiPatch.start();

                            neighbourFace[nNeighbourItems] = neighbourFaceI;

                            const vector curDelta =
                                neiPatch.neighbFaceCellCentres()[neighbourFaceI]
                                -
                                acceptorPos;

                            register const scalar w =
                                1.0/(SMALL + mag( curDelta ) );
                            neighbourWeights[nNeighbourItems] = w;
                            weightsum += w;

                            neighbourProc[nNeighbourItems]
                                = procID;
                        }
                    }
                }
            }
            primaryWeights[donorI] /= weightsum;

            if(nItems != 0 || nNeighbourItems != 0)
            {
                FatalErrorIn("void Foam::donorCellTetInterpolationMethod::update()")
                <<"Could not find all interpolation partners."
                <<abort(FatalError);
            }

            forAll(weights, itemI)
            {
                ownItems[nOwnInterpolationItems++].set
                (
                    donorI,
                    donors[itemI],
                    weights[itemI]/weightsum
                );
            }

            forAll(neighbourWeights, itemI)
            {
                const label procID = neighbourProc[itemI];
                const label neighbourItemI = nNeighbourItemsPerProc[procID]++;
                neighbourItems[procID][neighbourItemI].set
                (
                    donorI,
                    neighbourFace[itemI],
                    neighbourWeights[itemI]/weightsum
                );
                neighbourValueToFieldMap[procID][neighbourItemI]=donorI;
            }
        }
    }

    distribute(procToPatchMap);
}

