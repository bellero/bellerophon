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
    Foam::bellerophonNoPreconditioner

Description
    Null preconditioner for bellerophonLduMatrices.

SourceFiles
    bellerophonNoPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonPreconditioner_H
#define bellerophonPreconditioner_H

#include "bellerophonLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class bellerophonNoPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class bellerophonLduMatrix::preconditioner
{
    protected:

        // Protected data

            //- Reference to the base-solver this preconditioner is used with
            const bellerophonLduMatrix& matrix_;


    public:

        //- Find the preconditioner name (directly or from a sub-dictionary)
        static word getName(const dictionary&);

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                bellerophonMatrix,
                (
                    const bellerophonLduMatrix& mat,
                    const dictionary& solverControls
                ),
                (mat, solverControls)
            );

        // Constructors

            preconditioner
            (
                const bellerophonLduMatrix& mat
            )
            :
                matrix_(mat)
            {}


        // Selectors

            //- Return a new preconditioner
            static autoPtr<preconditioner> New
            (
                const bellerophonLduMatrix& mat,
                const dictionary& solverControls
            );


        //- Destructor
        virtual ~preconditioner()
        {}


        // Member functions

            //- Read and reset the preconditioner parameters
            //  from the given stream
            virtual void read(const dictionary&)
            {}

            //- Return wA the preconditioned form of residual rA
            virtual void precondition
            (
                scalarField& wA,
                const scalarField& rA,
                const direction cmpt=0
            ) const = 0;

            //- Return wT the transpose-matrix preconditioned form of
            //  residual rT.
            //  This is only required for preconditioning asymmetric matrices.
            virtual void preconditionT
            (
                scalarField& wT,
                const scalarField& rT,
                const direction cmpt=0
            ) const
            {
                notImplemented
                (
                    type() +"::preconditionT"
                    "(scalarField& wT, const scalarField& rT, "
                    "const direction cmpt)"
                );
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
