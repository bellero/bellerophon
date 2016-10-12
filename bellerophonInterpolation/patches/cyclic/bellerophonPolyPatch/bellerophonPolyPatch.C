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

#include "bellerophonPolyPatch.H"
#include "transformField.H"
#include "SubField.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaIntersect.H"
#include "ops.H"
#include <../../OpenFOAM/primitives/Tuple2/Tuple2.H>
#include "syncTools.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "faceSet.H"
#include "cellSet.H"
#include "meshTriangulation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, bellerophonPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, bellerophonPolyPatch, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::bellerophonPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initGeometry(pBufs);
}


void Foam::bellerophonPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    polyPatch::calcGeometry(pBufs);
}


void Foam::bellerophonPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::bellerophonPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
}


void Foam::bellerophonPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::bellerophonPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::bellerophonPolyPatch::bellerophonPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word &patchType,
    const transformType transform
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    donorZoneName_(""),
    oversetZoneName_(""),
    oversetZoneID_(-1)
{}


Foam::bellerophonPolyPatch::bellerophonPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word &patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    donorZoneName_(dict.lookup("donorZone")),
    oversetZoneName_
        (dict.lookupOrDefault<word>("oversetZone", "")),
    oversetZoneID_(-1)
{}


Foam::bellerophonPolyPatch::bellerophonPolyPatch
(
    const bellerophonPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    donorZoneName_(pp.donorZoneName_),
    oversetZoneName_(pp.oversetZoneName_),
    oversetZoneID_(pp.oversetZoneID_)
{}


Foam::bellerophonPolyPatch::bellerophonPolyPatch
(
    const bellerophonPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    donorZoneName_(pp.donorZoneName_),
    oversetZoneName_(pp.oversetZoneName_),
    oversetZoneID_(pp.oversetZoneID_)
{}


Foam::bellerophonPolyPatch::bellerophonPolyPatch
(
    const bellerophonPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    donorZoneName_(pp.donorZoneName_),
    oversetZoneName_(pp.oversetZoneName_),
    oversetZoneID_(pp.oversetZoneID_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bellerophonPolyPatch::~bellerophonPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::bellerophonPolyPatch::oversetZoneID() const
{
    if(oversetZoneID_ == -1)
    {
        if(oversetZoneName_ == "")
        {
            label zoneID = -1;
            if(this->faceCells().size() > 0)
            {
                // If no cells have been defined in the dictionary, check of first cell
                // of cells behind the interface

                const label firstCell = this->faceCells()[0];
                zoneID =
                    this->boundaryMesh().mesh().cellZones().whichZone(firstCell);
            }
            const scalar zoneIDCpy = zoneID;
            Pstream::gather(zoneID, maxOp<label>());
            Pstream::scatter(zoneID);
            if(zoneID < 0)
            {
                FatalErrorIn
                (
                    "Foam::labelList "
                    "Foam::bellerophonPolyPatch::oversetZones()  const"
                )
                << "Cannot find cell zone for patch "<< this->name() << "."
                << abort(FatalError);
            }
            else if (zoneIDCpy != -1 && zoneIDCpy != zoneID)
            {
                FatalErrorIn
                (
                    "Foam::labelList "
                    "Foam::bellerophonPolyPatch::oversetZones()  const"
                )
                << "Different zones next to patch "<< this->name() << "." << nl
                << "Please use oversetNames to specify overset zone names for this "
                << "patch explicitly."
                << abort(FatalError);
            }

            oversetZoneID_ = zoneID;

            oversetZoneName_ =
                this->boundaryMesh().mesh().cellZones().names()[zoneID];
        }
        else
        {
            oversetZoneID_ = this->boundaryMesh().mesh().cellZones()
                        .findZoneID(oversetZoneName_);

            if(oversetZoneID_ < 0)
            {
                FatalErrorIn
                (
                    "Foam::bellerophonPolyPatch::bellerophonPolyPatch "
                    "( "
                    "    const word& name, "
                    "    const dictionary& dict, "
                    "    const label index, "
                    "    const polyBoundaryMesh& bm, "
                    "    const word &patchType "
                    ")"
                )
                << "Cannot find cell zone "<<oversetZoneName_ << "."
                << abort(FatalError);
            }
        }
    }

    return oversetZoneID_;
}

Foam::autoPtr<Foam::triSurface> Foam::bellerophonPolyPatch::holeBoundary() const
{
    const polyMesh& mesh = this->boundaryMesh().mesh();

    const label nFaces = mesh.nFaces();
    const label nIntFaces = mesh.nInternalFaces();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    const labelList& pcs = this->faceCells();

    // Mark cells next to the patch and neighbours of these patches
    labelList marked(mesh.nCells(),0);
    {
        forAll(pcs, pcI)
        {
            marked[pcs[pcI]] = 1;
        }

        forAll(nei, faceI)
        {
            const label o = own[faceI];
            const label n = nei[faceI];
            if(marked[o] == 1 && marked[n] == 0)
            {
                marked[n] = 2;
            }
            else if(marked[o] == 0 && marked[n] == 1)
            {
                marked[o] = 2;
            }
        }
    }

    // Mark faces between marked and not marked cells
    boolList interfaceFace(nFaces,false);
    {
        for(label faceI = 0; faceI<nIntFaces; faceI++)
        {
            const label o = own[faceI];
            const label n = nei[faceI];
            if((marked[o] > 0) != (marked[n] > 0))
            {
                interfaceFace[faceI]=true;
            }
        }

        // Sync if parallel
        if(Pstream::parRun())
        {
            // This is a bit tricky as there is no notEqOp

            // Copy state to separate list to swap
            boolList neighbourInterfaceFace
                (nFaces-nIntFaces, false);

            const polyBoundaryMesh& patches = mesh.boundaryMesh();

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
                    const label start = patch.start() - nIntFaces;
                    for
                    (
                        label boundaryI = start;
                        boundaryI < start + patch.size();
                        boundaryI++
                    )
                    {
                        const label faceI = boundaryI + nIntFaces;
                        if(marked[own[faceI]] > 0)
                        {
                            interfaceFace[faceI] = true;
                            neighbourInterfaceFace[boundaryI] = true;
                        }
                    }
                }
            }

            // Swap state
            syncTools::swapBoundaryFaceList(mesh, neighbourInterfaceFace);

            // If only side marks boundary face it is part of the interface,
            // otherwise not
            for(label faceI = 0; faceI < nFaces-nIntFaces; faceI++)
            {
                interfaceFace[faceI+nIntFaces] =
                (
                    interfaceFace[faceI+nIntFaces]
                    !=
                    neighbourInterfaceFace[faceI]
                );
            }
        }

    }

    // Count interface triangles
    label nTris = 0;

    // Count face centre triangles
    label nCentres = 0;

    // Faces of the mesh
    const faceList& fs = mesh.faces();

    // Vertices of the mesh
    const pointField& ps = mesh.points();

    // Vertices of the mesh
    const pointField& fcs = mesh.faceCentres();

    // Count supporting points of triSurfaces
    label nMeshPoints = 0;

    // Mark triSurface point indices
    labelList triPointIndex(mesh.nPoints(), -1);

    forAll(interfaceFace, faceI)
    {
        if(interfaceFace[faceI])
        {
            const face& f = fs[faceI];
            nTris += f.size();
            nCentres ++;
            forAll(f, i)
            {
                if(triPointIndex[f[i]] == -1)
                {
                    triPointIndex[f[i]] = nMeshPoints;
                    nMeshPoints++;
                }
            }
        }
    }

    // Triangles
    List<labelledTri> tris(nTris);

    // Number of vertices supporting the triangulated surface
    const label nLocalTotalPoints = nMeshPoints + nCentres;

    // Vertices supporting the triangulated surface
    pointField points(nLocalTotalPoints);

    // Fill vertices
    forAll(triPointIndex, pointI)
    {
        if(triPointIndex[pointI] > -1)
        {
            points[triPointIndex[pointI]] = ps[pointI];
        }
    }

    // Count triangles
    label triI = 0;

    // Offset for centre vertex of face triangulation
    label centreI = nMeshPoints;

    // Fill triangles
    forAll(interfaceFace, faceI)
    {
        if(interfaceFace[faceI])
        {
            const face& f = fs[faceI];
            points[centreI] = fcs[faceI];
            label prev = f[f.size()-1];
            forAll(f, i)
            {
                label cur = f[i];
                // Cells outside of the hole boundary are marked, so flip
                // triangle if owner is marked
                if(marked[own[faceI]] > 0)
                {
                    tris[triI] =
                        labelledTri
                        (
                            triPointIndex[cur],
                            triPointIndex[prev],
                            centreI,
                            0
                        );
                }
                else
                {
                    tris[triI] =
                        labelledTri
                        (
                            triPointIndex[prev],
                            triPointIndex[cur],
                            centreI,
                            0
                        );
                }
                prev = cur;
                triI++;
            }
            centreI++;
        }
    }

    if(Pstream::parRun())
    {
        // Collect global number of points
        List< Tuple2<label, label> > pointsPerProc(Pstream::nProcs());

        pointsPerProc[Pstream::myProcNo()] =
            Tuple2<label, label>(Pstream::myProcNo(),nLocalTotalPoints);

        Pstream::gatherList(pointsPerProc);
        Pstream::scatterList(pointsPerProc);

        // Offset for point numbering
        label myStart = 0;

        // Global number of points
        label nTotal = 0;

        forAll(pointsPerProc, I)
        {
            nTotal += pointsPerProc[I].second();
            if(pointsPerProc[I].first() < Pstream::myProcNo())
            {
                myStart += pointsPerProc[I].second();
            }
        }

        // Resize position points and renumber triangles
        points.setSize(nTotal, vector::zero);
        if(myStart > 0)
        {
            // Move points in list if necessary
            for(label pointI = nLocalTotalPoints-1; pointI >= 0; pointI--)
            {
                points[pointI+myStart] = points[pointI];
                points[pointI] = vector::zero;
            }

            // Update tri labels
            forAll(tris, triI)
            {
                labelledTri& curTri = tris[triI];
                curTri[0]+=myStart;
                curTri[1]+=myStart;
                curTri[2]+=myStart;
            }
        }

        // Exchange points
        Pstream::listCombineGather(points,plusEqOp<vector>());
        Pstream::listCombineScatter(points);

        // Collect triangles
        List< List<labelledTri> > trisPerProc(Pstream::nProcs());
        trisPerProc[Pstream::myProcNo()] = tris;

        Pstream::gatherList(trisPerProc);
        Pstream::scatterList(trisPerProc);

        // Count global triangles
        label triI = 0;
        forAll(trisPerProc, procI)
        {
            triI += trisPerProc[procI].size();
        }

        // Combine triangles
        tris.setSize(triI);
        triI = 0;
        forAll(trisPerProc, procI)
        {
            const List<labelledTri> curTris = trisPerProc[procI];
            forAll(curTris, curTriI)
            {
                tris[triI] = curTris[curTriI];
                triI++;
            }
        }
    }

    autoPtr<triSurface> resultPtr(new triSurface(tris, points));

//     resultPtr().write(this->name()+"_hole"+Foam::name(Pstream::myProcNo())+"_"+mesh.time().timeName()+".stl");

    return resultPtr;
}


void Foam::bellerophonPolyPatch::transformPosition(pointField& l) const
{

    if (separated())
        {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract

        const vectorField& s = separation();
        if (s.size() == 1)
        {
            forAll(l, i)
            {
                l[i] -= s[0];
            }
        }
        else
        {
            l -= s;
        }
    }
}


void Foam::bellerophonPolyPatch::transformPosition
(
    point& l,
    const label faceI
) const
{}

void Foam::bellerophonPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{}


void Foam::bellerophonPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::bellerophonPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    // do nothing
    return false;
}


void Foam::bellerophonPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    os.writeKeyword("donorZone")<<donorZoneName_
            << token::END_STATEMENT << nl;

    if(oversetZoneName_ != "")
    {
        os.writeKeyword("oversetZone") << oversetZoneName_
            << token::END_STATEMENT << nl;
    }
}

// ************************************************************************* //
