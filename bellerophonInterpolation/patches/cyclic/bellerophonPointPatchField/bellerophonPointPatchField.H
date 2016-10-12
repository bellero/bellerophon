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
    Foam::bellerophonPointPatchField

Description
    Cyclic bellerophon front and back plane patch field

SourceFiles
    bellerophonPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonPointPatchField_H
#define bellerophonPointPatchField_H

#include "coupledPointPatchField.H"
#include "bellerophonPointPatch.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class bellerophonPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bellerophonPointPatchField
:
    public coupledPointPatchField<Type>
{
    // Private data

        //- Local reference cast into the bellerophon patch
        const bellerophonPointPatch& bellerophonPatch_;

        //- Owner side patch interpolation pointer
        mutable autoPtr<PrimitivePatchInterpolation<primitivePatch> > ppiPtr_;


    // Private Member Functions

        //- Owner side patch interpolation
        const PrimitivePatchInterpolation<primitivePatch>& ppi() const
        {
            Info<< "bellerophonPointPatchField::ppi()"
                << " at bellerophonPointPatchField.H:72" << nl << endl;

            if (!ppiPtr_.valid())
            {
                ppiPtr_.reset
                (
                    new PrimitivePatchInterpolation<primitivePatch>
                    (
                        bellerophonPatch_.bellerophonPatch()
                    )
                );
            }

            return ppiPtr_();
        }


public:

    //- Runtime type information
    TypeName(bellerophonPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        bellerophonPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bellerophonPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        bellerophonPointPatchField
        (
            const bellerophonPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type> > clone() const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new bellerophonPointPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        bellerophonPointPatchField
        (
            const bellerophonPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type> > clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new bellerophonPointPatchField<Type>
                (
                    *this, iF
                )
            );
        }


    // Member functions

        // Constraint handling

            //- Return the constraint type this coupledPointPatchField implements
            virtual const word& constraintType() const
            {
                Info<<"bellerophonPointPatchField::constraintType()"<<nl<<endl;

                return bellerophonPointPatch::typeName;
            }


/*
 *

        // Cyclic coupled interface functions

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return
                   !(
                        bellerophonPatch_.parallel()
                     || pTraits<Type>::rank == 0
                    );
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                Info<< "bellerophonPointPatchField::forwardT()"
                    << " at bellerophonPointPatchField.H:215" << nl << endl;

                return bellerophonPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                Info<< "bellerophonPointPatchField::reverseT()"
                    << " at bellerophonPointPatchField.H:224" << nl << endl;

                return bellerophonPatch_.reverseT();
            }
 */

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            )
            {}

            //- Complete swap of patch point values and add to local values
            virtual void swapAddSeparated
            (
                const Pstream::commsTypes commsType,
                Field<Type>&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "bellerophonPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
