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

#include "bellerophonDiagonalPreconditioner.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonDiagonalPreconditioner, 0);

    bellerophonLduMatrix::preconditioner::
        addbellerophonMatrixConstructorToTable<bellerophonDiagonalPreconditioner>
        addbellerophonDiagonalPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bellerophonDiagonalPreconditioner::bellerophonDiagonalPreconditioner
(
    const bellerophonLduMatrix& mat,
    const dictionary&
)
:
    bellerophonLduMatrix::preconditioner(mat),
    rD(mat.diag().size())
{
    scalar* __restrict__ rDPtr = rD.begin();
    const scalar* __restrict__ DPtr = matrix_.diag().begin();

    register label nCells = rD.size();

    // Generate reciprocal diagonal
    for (register label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/DPtr[cell];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bellerophonDiagonalPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction
) const
{
    scalar* __restrict__ wAPtr = wA.begin();
    const scalar* __restrict__ rAPtr = rA.begin();
    const scalar* __restrict__ rDPtr = rD.begin();

    register label nCells = wA.size();

    for (register label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }
}


// ************************************************************************* //
