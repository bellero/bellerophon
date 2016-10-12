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

#ifndef bellerophonLduMatrixI_H
#define bellerophonLduMatrixI_H


#include "bellerophonLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalarField& Foam::bellerophonLduMatrix::modLower()
{
    return *modLowerPtr_;
}

inline Foam::scalarField& Foam::bellerophonLduMatrix::modUpper()
{
    return *modUpperPtr_;
}

inline const Foam::scalarField& Foam::bellerophonLduMatrix::modLower() const
{
    return *modLowerPtr_;
}

inline const Foam::scalarField& Foam::bellerophonLduMatrix::modUpper() const
{
    return *modUpperPtr_;
}

inline const Foam::labelListList& Foam::bellerophonLduMatrix::donorCol() const
{
    return donorCols_;
};

inline const Foam::labelListList& Foam::bellerophonLduMatrix::acceptorRow() const
{
    return acceptorRows_;
};

#endif // bellerophonLduMatrixI_H