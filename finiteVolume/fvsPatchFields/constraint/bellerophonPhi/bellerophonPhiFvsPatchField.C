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

#include "bellerophonPhiFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::bellerophonPhiFvsPatchField<Type>::bellerophonPhiFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    bellerophonPatch_(refCast<const bellerophonFvPatch>(p))
{}


template<class Type>
Foam::bellerophonPhiFvsPatchField<Type>::bellerophonPhiFvsPatchField
(
    const bellerophonPhiFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    bellerophonPatch_(refCast<const bellerophonFvPatch>(p))
{
    if (!isA<bellerophonFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "bellerophonPhiFvsPatchField<Type>::bellerophonPhiFvsPatchField\n"
            "("
                "const bellerophonPhiFvsPatchField<Type>&, "
                "const fvPatch&, "
                "const DimensionedField<Type, surfaceMesh>&"
            ")"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::bellerophonPhiFvsPatchField<Type>::bellerophonPhiFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    bellerophonPatch_(refCast<const bellerophonFvPatch>(p))
{
    if (!isA<bellerophonFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "bellerophonPhiFvsPatchField<Type>::bellerophonPhiFvsPatchField"
            "("
                "const fvPatch&, "
                "const Field<Type>&, "
                "const dictionary&"
            ")",
            dict
        )   << "patch " << this->patch().index() << " not bellerophon type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::bellerophonPhiFvsPatchField<Type>::bellerophonPhiFvsPatchField
(
    const bellerophonPhiFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    bellerophonPatch_(ptf.bellerophonPatch_)
{}


template<class Type>
Foam::bellerophonPhiFvsPatchField<Type>::bellerophonPhiFvsPatchField
(
    const bellerophonPhiFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    bellerophonPatch_(ptf.bellerophonPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
namespace Foam
{

template<class Type>
void bellerophonPhiFvsPatchField<Type>::correctSum()
{
//     // Integrate over faces and sum up face are
//     // TODO: make this demand driven, so it is saved and updated on demand
// 
//     const scalarField areas = mag(bellerophonPatch_.Sf());
//     scalar area = 0.0;
//     Type sum = pTraits<Type>::zero;
// 
//     forAll(areas, areaI)
//     {
//         area+=areas[areaI];
//         sum+=areas[areaI]*(this->operator[](areaI));
//     }
// 
//     // Calculate specific continuity error
//     sum/=area;
// 
//     Info<<"Correcting sum, specific error is "<<areasum<<nl<<endl;
// 
//     // Correct values
//     forAll(areas, areaI)
//     {
//         this->operator[](areaI)-=sum*areas[areaI];
//     }
}

template<>
void bellerophonPhiFvsPatchField<scalar>::correctSum()
{
    // Integrate over faces and sum up face are
    // TODO: make this demand driven, so it is saved and updated on demand

    Info<<"Correcting field "<<*this<<endl;

    const scalarField areas = mag(bellerophonPatch_.Sf());
    scalar area = 0.0;
    scalar areasum = 0.0;

    forAll(areas, areaI)
    {
        area+=areas[areaI];
        areasum+=(this->operator[](areaI));
    }

    // Calculate specific continuity error
    areasum/=area;

//     Info<<"Correcting sum, specific error is "<<areasum<<endl;
//     scalar check = 0.0;
    // Correct values
    forAll(areas, areaI)
    {
        this->operator[](areaI)-=areasum*areas[areaI];
//         check+=this->operator[](areaI);
    }

//     Info<<"    after correction: "<<check<<nl<<endl;
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator=
(
    const UList<Type>& rhs
)
{
    fvsPatchField<Type>::operator=(rhs);

    Info<<"operator="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator=
(
    const fvsPatchField<Type>& rhs
)
{
    fvsPatchField<Type>::operator=(rhs);

    Info<<"operator="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator+=
(
    const fvsPatchField<Type>& rhs
)
{
    fvsPatchField<Type>::operator+=(rhs);

    Info<<"operator+="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator-=
(
    const fvsPatchField<Type>& rhs
)
{
    fvsPatchField<Type>::operator-=(rhs);

    Info<<"operator-="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator*=
(
    const fvsPatchField<scalar>& rhs
)
{
    fvsPatchField<Type>::operator*=(rhs);

    Info<<"operator*="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator/=
(
    const fvsPatchField<scalar>& rhs
)
{
    fvsPatchField<Type>::operator/=(rhs);

    Info<<"operator/="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator+=
(
    const Field<Type>& rhs
)
{
    fvsPatchField<Type>::operator+=(rhs);

    Info<<"operator+="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator-=
(
    const Field<Type>& rhs
)
{
    fvsPatchField<Type>::operator-=(rhs);

    Info<<"operator-="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator*=
(
    const Field<scalar>& rhs
)
{
    fvsPatchField<Type>::operator*=(rhs);

    Info<<"operator*="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator/=
(
    const Field<scalar>& rhs
)
{
    fvsPatchField<Type>::operator/=(rhs);

    Info<<"operator/="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator=
(
    const Type& rhs
)
{
    fvsPatchField<Type>::operator=(rhs);

    Info<<"operator="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator+=
(
    const Type& rhs
)
{
    fvsPatchField<Type>::operator+=(rhs);

    Info<<"operator+="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator-=
(
    const Type& rhs
)
{
    fvsPatchField<Type>::operator-=(rhs);

    Info<<"operator-="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator*=
(
    const scalar rhs
)
{
    fvsPatchField<Type>::operator*=(rhs);

    Info<<"operator*="<<endl;
    correctSum();
}

template<class Type>
void bellerophonPhiFvsPatchField<Type>::operator/=
(
    const scalar rhs
)
{
    fvsPatchField<Type>::operator/=(rhs);

    Info<<"operator/="<<endl;
    correctSum();
}

}

// ************************************************************************* //
