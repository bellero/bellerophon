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
#include "staticOversetFvMesh.H"

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

#include "bellerophonInterpolation.H"
#include "bellerophonInterface.H"
#include "bellerophonPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(staticOversetFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, staticOversetFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::staticOversetFvMesh::staticOversetFvMesh(const IOobject& io)
:
    solidOversetMotionFvMesh(io)
{
    this->moving(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::staticOversetFvMesh::~staticOversetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::staticOversetFvMesh::update()
{
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
                    "bool Foam::staticOversetFvMesh::update()"
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
                    const bool   flip  = interfaceFlipMap[faceI];
                    const label  cellI = flip ? nei[f] : own[f];
                    faceSums[cellI] += intMagSf[faceI];
                }

                forAll(interfaceFaces, faceI)
                {
                    const label  f     = interfaceFaces[faceI];
                    const bool   flip  = interfaceFlipMap[faceI];
                    const label  cellI = flip ? nei[f] : own[f];

                    const vector n     = intSf[f] / intMagSf[f];
                    const scalar vol   = cellVols[cellI];
                    const scalar area   =
                        flip ? faceSums[cellI] : -faceSums[cellI];
                    const vector factor = vol * n / area;

                    intUf[f] += intDivPhi[cellI] * factor;
                    intUfOld[f] = intUf[f];
                }
            }
        }

        // Now correct boundary conditions of velocity field
        const_cast<volVectorField&>(U).correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningIn("staticOversetFvMesh::update()")
            << "Could not find volVectorField " << UName_
            << ". Therefore not going to updating boundary conditions."
            << endl;
    }

    return false;
}

bool Foam::staticOversetFvMesh::writeObject
(
    Foam::IOstream::streamFormat fmt,
 Foam::IOstream::versionNumber ver,
 Foam::IOstream::compressionType cmp
) const
{
    bool writtenLiveCells = writeLiveCells(this->time().timeName());

    return writtenLiveCells && fvMesh::writeObject(fmt, ver, cmp);
}

// ************************************************************************* //
