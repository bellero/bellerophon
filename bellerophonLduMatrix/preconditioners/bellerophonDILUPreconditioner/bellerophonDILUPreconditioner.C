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

#include "bellerophonDILUPreconditioner.H"

#include "bellerophonPBiCG.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonDILUPreconditioner, 0);

    bellerophonLduMatrix::preconditioner::
        addbellerophonMatrixConstructorToTable<bellerophonDILUPreconditioner>
        addbellerophonDILUPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bellerophonDILUPreconditioner::bellerophonDILUPreconditioner
(
    const bellerophonLduMatrix& mat,
    const dictionary&
)
:
    bellerophonLduMatrix::preconditioner(mat),
    rD_(mat.diag())
{
   calcReciprocalD(rD_, mat);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bellerophonDILUPreconditioner::calcReciprocalD
(
    scalarField& rD,
    const bellerophonLduMatrix& mat
)
{
    scalar* __restrict__ rDPtr = rD.begin();

    const label* const __restrict__ uPtr = mat.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = mat.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = mat.modUpper().begin();
    const scalar* const __restrict__ lowerPtr = mat.modLower().begin();

    register label nFaces = mat.modUpper().size();
    for (register label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= upperPtr[face]*lowerPtr[face]/rDPtr[lPtr[face]];
    }


    // Calculate the reciprocal of the preconditioned diagonal
    register label nCells = rD.size();

    for (register label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }
}


void Foam::bellerophonDILUPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction
) const
{
    scalar* __restrict__ wAPtr = wA.begin();
    const scalar* __restrict__ rAPtr = rA.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();
    const label* const __restrict__ losortPtr =
        matrix_.lduAddr().losortAddr().begin();

    const scalar* const __restrict__ upperPtr =
        matrix_.modUpper().begin();
    const scalar* const __restrict__ lowerPtr =
        matrix_.modLower().begin();

    register label nCells = wA.size();
    register label nFaces = matrix_.modUpper().size();
    register label nFacesM1 = nFaces - 1;

    for (register label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }


    register label sface;

    for (register label face=0; face<nFaces; face++)
    {
        sface = losortPtr[face];
        wAPtr[uPtr[sface]] -=
            rDPtr[uPtr[sface]]*lowerPtr[sface]*wAPtr[lPtr[sface]];
    }

    for (register label face=nFacesM1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -=
            rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
    }
}


void Foam::bellerophonDILUPreconditioner::preconditionT
(
    scalarField& wT,
    const scalarField& rT,
    const direction
) const
{
    scalar* __restrict__ wTPtr = wT.begin();
    const scalar* __restrict__ rTPtr = rT.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();
    const label* const __restrict__ losortPtr =
        matrix_.lduAddr().losortAddr().begin();

    const scalar* const __restrict__ upperPtr =
        matrix_.modUpper().begin();
    const scalar* const __restrict__ lowerPtr =
        matrix_.modLower().begin();

    register label nCells = wT.size();
    register label nFaces = matrix_.modUpper().size();
    register label nFacesM1 = nFaces - 1;

    for (register label cell=0; cell<nCells; cell++)
    {
        wTPtr[cell] = rDPtr[cell]*rTPtr[cell];
    }

    for (register label face=0; face<nFaces; face++)
    {
        wTPtr[uPtr[face]] -=
            rDPtr[uPtr[face]]*upperPtr[face]*wTPtr[lPtr[face]];
    }


    register label sface;

    for (register label face=nFacesM1; face>=0; face--)
    {
        sface = losortPtr[face];
        wTPtr[lPtr[sface]] -=
            rDPtr[lPtr[sface]]*lowerPtr[sface]*wTPtr[uPtr[sface]];
    }
}


// ************************************************************************* //
