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

#include "bellerophonPointPatchField.H"
#include "Swap.H"
#include "transformField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::bellerophonPointPatchField<Type>::bellerophonPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    bellerophonPatch_(refCast<const bellerophonPointPatch>(p)),
    ppiPtr_(NULL)
{}


template<class Type>
Foam::bellerophonPointPatchField<Type>::bellerophonPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    bellerophonPatch_(refCast<const bellerophonPointPatch>(p)),
    ppiPtr_(NULL)
{
    if (!isType<bellerophonPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "bellerophonPointPatchField<Type>::bellerophonPointPatchField\n"
            "(\n"
            "    const pointPatch&,\n"
            "    const DimensionedField<Type, pointMesh>&,\n"
            "    const dictionary&\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not bellerophon type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::bellerophonPointPatchField<Type>::bellerophonPointPatchField
(
    const bellerophonPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(/*ptf, */p, iF/*, mapper*/),
    bellerophonPatch_(refCast<const bellerophonPointPatch>(p)),
    ppiPtr_(NULL)
{
    if (!isType<bellerophonPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "bellerophonPointPatchField<Type>::bellerophonPointPatchField\n"
            "(\n"
            "    const bellerophonPointPatchField<Type>&,\n"
            "    const pointPatch&,\n"
            "    const DimensionedField<Type, pointMesh>&"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::bellerophonPointPatchField<Type>::bellerophonPointPatchField
(
    const bellerophonPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    bellerophonPatch_(ptf.bellerophonPatch_),
    ppiPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



template<class Type>
void Foam::bellerophonPointPatchField<Type>::swapAddSeparated
(
    const Pstream::commsTypes,
    Field<Type>& pField
) const
{
/*
    // We inplace modify pField. To prevent the other side (which gets
    // evaluated at a later date) using already changed values we do
    // all swaps on the side that gets evaluated first.

    // Get neighbouring pointPatch
    const bellerophonPointPatch& nbrPatch = bellerophonPatch_.neighbPatch();

    // Get neighbouring coupledPointPatchField
    const GeometricField<Type, coupledPointPatchField, pointMesh>& fld =
        refCast<const GeometricField<Type, coupledPointPatchField, pointMesh> >
        (
            this->dimensionedInternalField()
        );

    const bellerophonPointPatchField<Type>& nbr =
        refCast<const bellerophonPointPatchField<Type> >
        (
            fld.boundaryField()[nbrPatch.index()]
        );


    Field<Type> ptFld(this->patchInternalField(pField));
    Field<Type> nbrPtFld(nbr.patchInternalField(pField));


//     if (doTransform())
//     {
//         const tensor& forwardT = this->forwardT()[0];
//         const tensor& reverseT = this->reverseT()[0];
//
//         transform(ptFld, reverseT, ptFld);
//         transform(nbrPtFld, forwardT, nbrPtFld);
//     }

    // convert point field to face field, AMI interpolate, then
    // face back to point
    {
        // add neighbour side contribution to owner
        Field<Type> nbrFcFld(nbrPpi().pointToFaceInterpolate(nbrPtFld));

        // interpolate to owner
        nbrFcFld = bellerophonPatch_.bellerophonPatch().interpolate(nbrFcFld);

        // add to internal field
        this->addToInternalField
        (
            pField,
            ppi().faceToPointInterpolate(nbrFcFld)()
        );
    }
*/

}

// ************************************************************************* //
