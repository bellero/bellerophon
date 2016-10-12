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
    Foam::bellerophonGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface field.

SourceFiles
    bellerophonGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonGAMGInterfaceField_H
#define bellerophonGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "bellerophonGAMGInterface.H"
#include "bellerophonInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class bellerophonGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class bellerophonGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public bellerophonInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const bellerophonGAMGInterface& bellerophonInterface_;

/*
        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;
 */

    // Private Member Functions

        //- Disallow default bitwise copy construct
        bellerophonGAMGInterfaceField(const bellerophonGAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const bellerophonGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("bellerophon");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        bellerophonGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );


    //- Destructor
    virtual ~bellerophonGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                Info<< "bellerophonGAMGInterfaceField::size()"
                    << " at bellerophonGAMGInterfaceField.C:76" << nl << endl;

                return bellerophonInterface_.size();
            }


        // Interface matrix update

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

/*
 * stfu with you transformation

        //- Cyclic interface functions

            //- Does the interface field perform the transfromation
            virtual bool doTransform() const
            {
                Info<< "bellerophonGAMGInterfaceField::doTransform(...)"
                    << " at bellerophonGAMGInterfaceField.C:128" << nl << endl;

                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                Info<< "bellerophonGAMGInterfaceField::forwardT()"
                    << " at bellerophonGAMGInterfaceField.C:137" << nl << endl;

                return bellerophonInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                Info<< "bellerophonGAMGInterfaceField::reverseT()"
                    << " at bellerophonGAMGInterfaceField.C:146" << nl << endl;

                return bellerophonInterface_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                Info<< "bellerophonGAMGInterfaceField::rank()"
                    << " at bellerophonGAMGInterfaceField.C:155" << nl << endl;

                return rank_;
            }
*/

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
