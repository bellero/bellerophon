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

Class
    Foam::bellerophonPhiFvsPatchField

Description
    Foam::bellerophonPhiFvsPatchField

SourceFiles
    bellerophonPhiFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonPhiFvsPatchField_H
#define bellerophonPhiFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "bellerophonFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bellerophonPhiFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bellerophonPhiFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const bellerophonFvPatch& bellerophonPatch_;

public:

    //- Runtime type information
    TypeName("bellerophonPhi");


    // Constructors

        //- Construct from patch and internal field
        bellerophonPhiFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bellerophonPhiFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given bellerophonPhiFvsPatchField onto a new patch
        bellerophonPhiFvsPatchField
        (
            const bellerophonPhiFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        bellerophonPhiFvsPatchField
        (
            const bellerophonPhiFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {

            return tmp<fvsPatchField<Type> >
            (
                new bellerophonPhiFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        bellerophonPhiFvsPatchField
        (
            const bellerophonPhiFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {

            return tmp<fvsPatchField<Type> >
            (
                new bellerophonPhiFvsPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        //- Correct field to make sum of all elements zero
        //  At the moment, this is only implemented for scalar fvsPatchFields
        //  because as far as I can see, we use this only for the face flux
        //  If not, we will have to do something else...
        //  TODO:
        //  make the patch call this function every time it recalculates
        //  the neighbour field or so.
        //  For now, this is quick and very dirty but perhaps works.
        void correctSum();


        // Member operators

        virtual void operator=(const UList<Type>&);
        virtual void operator=(const fvsPatchField<Type>&);
        virtual void operator+=(const fvsPatchField<Type>&);
        virtual void operator-=(const fvsPatchField<Type>&);
        virtual void operator*=(const fvsPatchField<scalar>&);
        virtual void operator/=(const fvsPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "bellerophonPhiFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
