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

#include "bellerophonPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(bellerophonLduMatrix::preconditioner, bellerophonMatrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::bellerophonLduMatrix::preconditioner::getName
(
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


Foam::autoPtr<Foam::bellerophonLduMatrix::preconditioner>
Foam::bellerophonLduMatrix::preconditioner::New
(
    const bellerophonLduMatrix& mat,
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    bellerophonMatrixConstructorTable::iterator constructorIter =
        bellerophonMatrixConstructorTablePtr_->find(name);

    if (constructorIter == bellerophonMatrixConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "bellerophonLduMatrix::preconditioner::New"
            "(const solver&, const dictionary&)",
            controls
        )   << "Unknown asymmetric matrix preconditioner "
            << name << nl << nl
            << "Valid asymmetric matrix preconditioners :" << endl
            << bellerophonMatrixConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<bellerophonLduMatrix::preconditioner>
    (
        constructorIter()
        (
            mat,
            controls
        )
    );
}


// ************************************************************************* //
