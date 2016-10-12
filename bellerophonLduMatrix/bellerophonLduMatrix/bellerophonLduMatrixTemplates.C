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
    bellerophonLduMatrix member H operations.
    Seems like they are never called, so only overwrite them and throw
    an error of the get called.

\*---------------------------------------------------------------------------*/

#include "bellerophonLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::bellerophonLduMatrix::H(const Field<Type>& psi) const
{
    tmp<Field<Type> > tHpsi
    (
        new Field<Type>(lduAddr().size(), pTraits<Type>::zero)
    );

    FatalErrorIn
        ("Foam::tmp<Foam::Field<Type> >"
         "Foam::bellerophonLduMatrix::H(const Field<Type>& psi) const")
        << "not implemented"
        << abort(FatalError);

    return tHpsi;
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::bellerophonLduMatrix::H(const tmp<Field<Type> >& tpsi) const
{
    tmp<Field<Type> > tHpsi(H(tpsi()));
    tpsi.clear();
    return tHpsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::bellerophonLduMatrix::faceH(const Field<Type>& psi) const
{
    tmp<Field<Type> > tfaceHpsi(new Field<Type> (lduMatrix::lower().size()));

    FatalErrorIn
        ("Foam::tmp<Foam::Field<Type> >"
         "Foam::bellerophonLduMatrix::faceH(const Field<Type>& psi) const")
        << "not implemented"
        << abort(FatalError);

    return tfaceHpsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::bellerophonLduMatrix::faceH(const tmp<Field<Type> >& tpsi) const
{
    tmp<Field<Type> > tfaceHpsi(faceH(tpsi()));
    tpsi.clear();
    return tfaceHpsi;
}


// ************************************************************************* //
