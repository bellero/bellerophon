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

#include "bellerophonGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        bellerophonGAMGInterfaceField,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bellerophonGAMGInterfaceField::bellerophonGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    bellerophonInterface_(refCast<const bellerophonGAMGInterface>(GAMGCp))
//     doTransform_(false),
//     rank_(0)
{
    Info<< "bellerophonGAMGInterfaceField::bellerophonGAMGInterfaceField(...)"
        << " at bellerophonGAMGInterfaceField.C:46" << nl << endl;

    const bellerophonInterfaceField& p =
        refCast<const bellerophonInterfaceField>(fineInterface);

//     doTransform_ = p.doTransform();
//     rank_ = p.rank();
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::bellerophonGAMGInterfaceField::~bellerophonGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bellerophonGAMGInterfaceField::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    Info<< "bellerophonGAMGInterfaceField::updateInterfaceMatrix(...)"
        << " at bellerophonGAMGInterfaceField.C:76" << nl << endl;

//     // Get neighbouring field
//     scalarField pnf
//     (
//         bellerophonInterface_.neighbPatch().interfaceInternalField(psiInternal)
//     );
// 
//     // Transform according to the transformation tensors
//     transformCoupleField(pnf, cmpt);
// 
//     if (bellerophonInterface_.owner())
//     {
//         pnf = bellerophonInterface_.bellerophon().interpolateToSource(pnf);
//     }
//     else
//     {
//         pnf = bellerophonInterface_.neighbPatch().bellerophon().interpolateToTarget(pnf);
//     }
// 
//     const labelUList& faceCells = bellerophonInterface_.faceCells();
// 
//     forAll(faceCells, elemI)
//     {
//         result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
//     }
}


// ************************************************************************* //
