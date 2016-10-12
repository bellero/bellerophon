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
    Foam::bellerophonSmoother

Description
    Abstract base class for smoothers for bellerophonLduMatrices.

SourceFiles
    bellerophonSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonSmoother_H
#define bellerophonSmoother_H

#include "bellerophonLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class bellerophonSmoother Declaration
\*---------------------------------------------------------------------------*/

class bellerophonLduMatrix::smoother
{
protected:

    // Protected data

        word fieldName_;
        const bellerophonLduMatrix& matrix_;
        const FieldField<Field, scalar>& interfaceBouCoeffs_;
        const FieldField<Field, scalar>& interfaceIntCoeffs_;
        const lduInterfaceFieldPtrsList& interfaces_;


public:

    //- Find the smoother name (directly or from a sub-dictionary)
    static word getName(const dictionary&);

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            smoother,
            bellerophonMatrix,
            (
                const word& fieldName,
                const bellerophonLduMatrix& matrix,
                const FieldField<Field, scalar>& interfaceBouCoeffs,
                const FieldField<Field, scalar>& interfaceIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces
            ),
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces
            )
        );


    // Constructors

        smoother
        (
            const word& fieldName,
            const bellerophonLduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Selectors

        //- Return a new smoother
        static autoPtr<smoother> New
        (
            const word& fieldName,
            const bellerophonLduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~smoother()
    {}


    // Member functions

        // Access

            const word& fieldName() const
            {
                return fieldName_;
            }

            const bellerophonLduMatrix& matrix() const
            {
                return matrix_;
            }

                const FieldField<Field, scalar>& interfaceBouCoeffs() const
                {
                    return interfaceBouCoeffs_;
                }

                const FieldField<Field, scalar>& interfaceIntCoeffs() const
                {
                    return interfaceIntCoeffs_;
                }

                const lduInterfaceFieldPtrsList& interfaces() const
                {
                    return interfaces_;
                }


        //- Smooth the solution for a given number of sweeps
        virtual void smooth
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt,
            const label nSweeps
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //