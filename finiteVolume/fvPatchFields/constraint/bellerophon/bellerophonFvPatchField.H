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
    Foam::bellerophonFvPatchField

Description
    Foam::bellerophonFvPatchField

SourceFiles
    bellerophonFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonFvPatchField_H
#define bellerophonFvPatchField_H

#include "coupledFvPatchField.H"
#include "bellerophonInterfaceField.H"
#include "bellerophonFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bellerophonFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bellerophonFvPatchField
:
    virtual public bellerophonInterfaceField,
    public coupledFvPatchField<Type>
{
private:
    // Private data

        //- Local reference cast into the cyclic patch
        const bellerophonInterface& interface_;

        //- Force interpolation on assignment
        bool forceInterpolation_;

        //- Is the field bounded?
        bool bound_;

//         //- Pretend to fix the value?
//         //  This is needed when using moving meshes when the correcting weights
//         //  for the transient correction of fluxes is done:
//         //  fvc::ddtCorr(U,Uf) calling fvcDdtPhiCoeff
//         //  TODO one could implement a own ddtScheme, that reimplements
//         //  fvcDdtPhiCoeff
//         bool fixesValue_;

public:

    //- Runtime type information
    TypeName(bellerophonFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        bellerophonFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bellerophonFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given bellerophonFvPatchField onto a new patch
        bellerophonFvPatchField
        (
            const bellerophonFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        bellerophonFvPatchField(const bellerophonFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new bellerophonFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        bellerophonFvPatchField
        (
            const bellerophonFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new bellerophonFvPatchField<Type>(*this, iF)
            );
        }


// Member functions

    // Access

        virtual const bellerophonInterface& interface() const
        {
            return interface_;
        }

        virtual const word& iFieldName() const;

        virtual bool fixesValue() const { return false /*fixesValue_*/;};

    // Evaluation functions

        //- Return true if coupled. Note that the underlying patch
        //  is not coupled() - the points don't align.
        virtual bool coupled() const;

        //- Return neighbour coupled internal cell data
        virtual tmp<Field<Type> > patchNeighbourField() const;

        //- Return surface normal gradient
        virtual tmp<Field<Type> > snGrad() const;

        //- Return surface normal gradient
        virtual tmp<Field<Type> > snGrad(const scalarField&) const;

        //- Initialise the evaluation of the patch field.
        //  Interpolation needs to be done here, as cells values are read by
        //  the processor patches during initialisation
        virtual void initEvaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the
        // evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type> > valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        // evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type> > valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        // evaluation of the gradient of this patchField
        virtual tmp<Field<Type> > gradientInternalCoeffs
        (
            const scalarField& deltaCoeffs
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        // evaluation of the gradient of this patchField
        virtual tmp<Field<Type> > gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        // evaluation of the gradient of this patchField
        virtual tmp<Field<Type> > gradientBoundaryCoeffs
        (
            const scalarField& deltaCoeffs
        ) const;

        //- Return the matrix source coefficients corresponding to the
        // evaluation of the gradient of this patchField
        virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;

        //- Update result field based on interface functionality
        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        //- Update result field based on interface functionality
        virtual void updateInterfaceMatrix
        (
            Field<Type>& result,
            const Field<Type>& psiInternal,
            const scalarField& coeffs,
            const Pstream::commsTypes commsType
        ) const;

    // Write
        virtual void write(Ostream& os) const;


    // Overwrite member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

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

#endif

// ************************************************************************* //
