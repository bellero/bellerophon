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
    Foam::printMatrix

Description
    Preconditioned bi-conjugate gradient solver for asymmetric lduMatrices
    using a run-time selectable preconditiioner.

SourceFiles
    printMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef printMatrix_H
#define printMatrix_H

#include "lduMatrix.H"
#include "bellerophonLduMatrix.H"
#include "bellerophonInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class printMatrix Declaration
\*---------------------------------------------------------------------------*/

class printMatrix
:
    public lduMatrix::solver
{
    // Private data

//         //- List of bellerophon patches
//         UPtrList<const bellerophonInterfaceField> belleroInterfaces_;

        //- Corrected Matrix
        bellerophonLduMatrix bellerophonMatrix_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        printMatrix(const printMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const printMatrix&);

public:

    //- Runtime type information
    TypeName("printMatrix");

    // Constructors

        //- Construct from matrix components and solver data stream
        printMatrix
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~printMatrix()
    {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;

        const bellerophonLduMatrix& bellerophonMatrix() const
        {
            return bellerophonMatrix_;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
