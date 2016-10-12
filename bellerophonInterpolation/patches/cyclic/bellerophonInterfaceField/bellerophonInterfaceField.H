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
    Foam::bellerophonInterfaceField

Description
    Abstract base class for cyclic bellerophon interfaces

SourceFiles
    bellerophonInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonInterfaceField_H
#define bellerophonInterfaceField_H

#include "primitiveFieldsFwd.H"
#include "bellerophonLduMatrix.H"
#include "bellerophonInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class bellerophonInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class bellerophonInterfaceField
{
public:

    //- Runtime type information
    TypeName("bellerophonInterfaceField");


    // Constructors

        //- Construct null
        bellerophonInterfaceField();

    //- Destructor
    virtual ~bellerophonInterfaceField();


    // Member Functions

        virtual const bellerophonInterface& interface() const = 0;

        virtual void initEvaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        ) = 0;

        virtual const word& iFieldName() const = 0;

/*
 * one nice day transformation will be implemented
        // Access

            //- Is the transform required
            virtual bool doTransform() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;


        //- Transform given patch internal field
        void transformCoupleField
        (
            scalarField& psiInternal,
            const direction cmpt
        ) const;
 */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
