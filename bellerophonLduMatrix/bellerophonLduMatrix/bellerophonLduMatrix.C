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

#include "lduMatrix.H"
#include "IOstreams.H"
#include "bellerophonLduMatrix.H"
#include "bellerophonInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonLduMatrix, 1);
}

const Foam::scalar Foam::bellerophonLduMatrix::great_ = 1.0e+20;
const Foam::scalar Foam::bellerophonLduMatrix::small_ = 1.0e-20;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::bellerophonLduMatrix::bellerophonLduMatrix(const Foam::bellerophonLduMatrix& A)
:
    lduMatrix(dynamicCast<const lduMatrix>(A)),
    modUpperPtr_(NULL),
    modLowerPtr_(NULL),
    primaryDonorCells_(A.primaryDonorCells_),
    primaryDonorWeights_(A.primaryDonorWeights_),
    ownInterpolationItems_(A.ownInterpolationItems_),
    neighbourInterpolationItems_(A.neighbourInterpolationItems_),
    neighbourValueToFieldMap_(A.neighbourValueToFieldMap_),
    interpolatedPsi_(A.primaryDonorCells_.size()),
    donorCols_(A.donorCols_),
    acceptorRows_(A.acceptorRows_)
{
    modUpperPtr_ = new scalarField(A.upper());
    modLowerPtr_ = new scalarField(A.lower());

    resizeBufs();
}


Foam::bellerophonLduMatrix::bellerophonLduMatrix(const Foam::lduMatrix& A)
:
    lduMatrix(A),
    modUpperPtr_(NULL),
    modLowerPtr_(NULL),
    primaryDonorCells_(bellerophon::Interpolation().primaryDonorCells()),
    primaryDonorWeights_(bellerophon::Interpolation().primaryDonorWeights()),
    ownInterpolationItems_
        (bellerophon::Interpolation().ownInterpolationItems()),
    neighbourInterpolationItems_
        (bellerophon::Interpolation().neighbourInterpolationItems()),
    neighbourValueToFieldMap_(bellerophon::Interpolation().neighbourValueToFieldMap()),
    interpolatedPsi_(primaryDonorCells_.size()),
    donorCols_(bellerophon::Interpolation().donorCols()),
    acceptorRows_(bellerophon::Interpolation().acceptorRows())
{
    // not creating additional coeffs, since they are created from
    // bellerophonLduInterfaceFields by the solver

    modUpperPtr_ = new scalarField(A.upper());
    modLowerPtr_ = new scalarField(A.lower());

    resizeBufs();
}


Foam::bellerophonLduMatrix::~bellerophonLduMatrix()
{
    delete modUpperPtr_;
    delete modLowerPtr_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

void Foam::bellerophonLduMatrix::resizeBufs() const
{
    const label nProcs = Pstream::nProcs();
    acceptorRowBuf_.setSize(nProcs);
    donorColBuf_.setSize(nProcs);
    neighbourInterpolationBuf_.setSize(nProcs);
    neighbourValueBuf_.setSize(nProcs);

    for(label procI = 0; procI < nProcs; procI++ )
    {
        acceptorRowBuf_[procI].setSize(acceptorRows_[procI].size());
        donorColBuf_[procI].setSize(donorCols_[procI].size());
        neighbourInterpolationBuf_[procI].setSize
            (neighbourInterpolationItems_[procI].size());
        neighbourValueBuf_[procI].setSize
            (neighbourValueToFieldMap_[procI].size());
    }
}


Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const Foam::bellerophonLduMatrix& bellerophonLdum)
{
    os << dynamicCast<const lduMatrix>(bellerophonLdum);

    if (bellerophonLdum.modLowerPtr_)
    {
        os  << "Modified lower triangle = "
            << *bellerophonLdum.modLowerPtr_
            << endl << endl;
    }

    if (bellerophonLdum.modUpperPtr_)
    {
        os  << "Modified Upper Coeffs = "
            << *bellerophonLdum.modUpperPtr_
            << endl;
    }

    os  << "Donor columns = "
        << bellerophonLdum.donorCol()
        << nl
        << "Acceptor rows = "
        << bellerophonLdum.acceptorRow()
        << endl;

    os.check("Ostream& operator<<(Ostream&, const bellerophonLduMatrix&");

    return os;
}


// ************************************************************************* //
