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
    Foam::bellerophonFvsPatchField

Description
    Foam::bellerophonFvsPatchField

SourceFiles
    bellerophonFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonFvsPatchField_H
#define bellerophonFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "bellerophonFvPatch.H"
#include "surfaceMesh.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bellerophonFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bellerophonFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const bellerophonFvPatch& bellerophonPatch_;

        //- Shall the continuity be enforced on this patch?
        bool enforceContinuity_;

        //- Sum of fluxes to enforce continuity (positive means outflow)
        scalar massFlux_;

        //- Name of the field
        word fieldName_;

public:

    //- Runtime type information
    TypeName(bellerophonFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        bellerophonFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bellerophonFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given bellerophonFvsPatchField onto a new patch
        bellerophonFvsPatchField
        (
            const bellerophonFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        bellerophonFvsPatchField
        (
            const bellerophonFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {

            return tmp<fvsPatchField<Type> >
            (
                new bellerophonFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        bellerophonFvsPatchField
        (
            const bellerophonFvsPatchField<Type>&,
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
                new bellerophonFvsPatchField<Type>(*this, iF)
            );
        }

        virtual bool fixesValue() const { return true; };


    // Member functions

        //- Correct field to make sum of all elements zero
        //  At the moment, this is only implemented for scalar fvsPatchFields
        //  because as far as I can see, we use this only for the face flux
        //  If not, we will have to do something else...
        void correctSum();

        //- Use name of GeometricField to enforce continuity of get name from
        //  dimensionedInternalField
        void enforceContinuity();

        //- Enforce continuity with given surfaceScalarField
        void enforceContinuity(const GeometricField<Type, fvsPatchField, surfaceMesh>&);

        //- Assign uniform value from operator= to internal hole interface faces
        void assignHoleInterface(const Type& rhs);

        // Constant access to mass flux
        scalar massFlux() const { return massFlux_; };

        // Non-constant access to mass flux
        scalar& massFlux() { return massFlux_; };

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

    // Write operation
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

