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

#include "bellerophonPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        bellerophonPointPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::bellerophonPointPatch::initGeometry(PstreamBuffers&)
{}


void Foam::bellerophonPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::bellerophonPointPatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::bellerophonPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::bellerophonPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
//     Info<<"bellerophonPointPatch::initUpdateMesh(..)"<<nl<<endl;

    facePointPatch::initUpdateMesh(pBufs);
//    bellerophonPointPatch::initGeometry(pBufs);
}


void Foam::bellerophonPointPatch::updateMesh(PstreamBuffers& pBufs)
{
//     Info<<"bellerophonPointPatch::updateMesh(..)"<<nl<<endl;

    facePointPatch::updateMesh(pBufs);
//    bellerophonPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bellerophonPointPatch::bellerophonPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    bellerophonPolyPatch_(refCast<const bellerophonPolyPatch>(patch))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bellerophonPointPatch::~bellerophonPointPatch()
{}


// ************************************************************************* //
