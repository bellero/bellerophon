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

Description
    bellerophonLduMatrix member operations.

\*---------------------------------------------------------------------------*/

#include "bellerophonLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::bellerophonLduMatrix::sumDiag()
{
    // Per definition sumDiag is zero for acceptor cells

    lduMatrix::sumDiag();

    const labelList& acceptorRows = acceptorRows_[Pstream::myProcNo()];
    const label* const __restrict__ rowPtr = acceptorRows.begin();
    const label nRows = acceptorRows.size();

    scalar* __restrict__ diagPtr = diag().begin();

    for(register label rowI = 0; rowI < nRows; rowI++)
    {
        diagPtr[rowPtr[rowI]] = 0.0;
    }
}


void Foam::bellerophonLduMatrix::negSumDiag()
{
    // Cache old diag, cause it will be modified with old upper and lower!

    const labelList& acceptorRows = acceptorRows_[Pstream::myProcNo()];
    const label nRows = acceptorRows.size();

    scalar* __restrict__ tmpDiagPtr = new scalar[nRows];
    const label* const __restrict__ rowPtr = acceptorRows.begin();

    {
        const scalar* const __restrict__ diagPtr = diag().begin();

        for(register label rowI = 0; rowI < nRows; rowI++)
        {
            tmpDiagPtr[rowI] = diagPtr[rowPtr[rowI]];
        }
    }

    lduMatrix::negSumDiag();

    scalar* __restrict__ diagPtr = diag().begin();

    for(register label rowI = 0; rowI < nRows; rowI++)
    {
        diagPtr[rowPtr[rowI]] = 2.0 * tmpDiagPtr[rowI];
    }
}


void Foam::bellerophonLduMatrix::sumMagOffDiag
(
    scalarField& sumOff
) const
{
    // Cache old diag, cause it will be modified with old upper and lower!

    const labelList& acceptorRows = acceptorRows_[Pstream::myProcNo()];
    const label nRows = acceptorRows.size();

    scalar* __restrict__ tmpDiagPtr = new scalar[nRows];
    const label* const __restrict__ rowPtr = acceptorRows.begin();

    {
        const scalar* const __restrict__ diagPtr = diag().begin();

        for(register label rowI = 0; rowI < nRows; rowI++)
        {
            tmpDiagPtr[rowI] = diagPtr[rowPtr[rowI]];
        }
    }

    lduMatrix::sumMagOffDiag(sumOff);

    scalar* __restrict__ diagPtr =
        const_cast<scalarField&>(lduMatrix::diag()).begin();

    for(register label rowI = 0; rowI < nRows; rowI++)
    {
        diagPtr[rowPtr[rowI]] = 2.0 * tmpDiagPtr[rowI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::bellerophonLduMatrix::operator=(const bellerophonLduMatrix& A)
{
    FatalErrorIn
    (
        "void Foam::bellerophonLduMatrix::operator="
        "( const bellerophonLduMatrix& )"
    )
    << "Arithmetric operations not defined for bellerophonLduMatrices."
    << abort(FatalError);
}


void Foam::bellerophonLduMatrix::negate()
{
    lduMatrix::negate();
}


void Foam::bellerophonLduMatrix::operator+=(const bellerophonLduMatrix& A)
{
    FatalErrorIn
    (
        "void Foam::bellerophonLduMatrix::operator="
        "( const bellerophonLduMatrix& )"
    )
    << "Arithmetric operations not defined for bellerophonLduMatrices."
    << abort(FatalError);
}


void Foam::bellerophonLduMatrix::operator-=(const bellerophonLduMatrix& A)
{
    FatalErrorIn
    (
        "void Foam::bellerophonLduMatrix::operator="
        "( const bellerophonLduMatrix& )"
    )
    << "Arithmetric operations not defined for bellerophonLduMatrices."
    << abort(FatalError);
}


void Foam::bellerophonLduMatrix::operator*=(const scalarField& sf)
{
    FatalErrorIn
    (
        "void Foam::bellerophonLduMatrix::operator="
        "( const bellerophonLduMatrix& )"
    )
    << "Arithmetric operations not defined for bellerophonLduMatrices."
    << abort(FatalError);
}


void Foam::bellerophonLduMatrix::operator*=(scalar s)
{
    FatalErrorIn
    (
        "void Foam::bellerophonLduMatrix::operator="
        "( const bellerophonLduMatrix& )"
    )
    << "Arithmetric operations not defined for bellerophonLduMatrices."
    << abort(FatalError);
}


Foam::tmp<Foam::scalarField > Foam::bellerophonLduMatrix::H1() const
{
    FatalErrorIn
    (
        "void Foam::bellerophonLduMatrix::operator="
        "( const bellerophonLduMatrix& )"
    )
    << "Arithmetric operations not defined for bellerophonLduMatrices."
    << abort(FatalError);

    tmp<scalarField > tH1
    (
        new scalarField(1, 0.0)
    );

    return tH1;
}


// ************************************************************************* //
