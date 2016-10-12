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
    Foam::bellerophonDiagonalPreconditioner

Description
    Diagonal preconditioner for both symmetric and asymmetric matrices.

    The reciprocal of the diagonal is calculated and stored for reuse
    because on most systems '*' is faster than '/'.

SourceFiles
    bellerophonDiagonalPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonDiagonalPreconditioner_H
#define bellerophonDiagonalPreconditioner_H

#include "bellerophonPreconditioner.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class bellerophonDiagonalPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class bellerophonDiagonalPreconditioner
:
    public bellerophonLduMatrix::preconditioner
{
    // Private data

        //- The reciprocal diagonal
        scalarField rD;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        bellerophonDiagonalPreconditioner(const bellerophonDiagonalPreconditioner&);

        //- Disallow default bitwise assignment
        void operator=(const bellerophonDiagonalPreconditioner&);


public:

    //- Runtime type information
    TypeName("diagonal");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        bellerophonDiagonalPreconditioner
        (
            const bellerophonLduMatrix& mat,
            const dictionary& solverControlsUnused
        );


    //- Destructor
    virtual ~bellerophonDiagonalPreconditioner()
    {}


    // Member Functions

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt=0
        ) const;

        //- Return wT the transpose-matrix preconditioned form of residual rT.
        virtual void preconditionT
        (
            scalarField& wT,
            const scalarField& rT,
            const direction cmpt=0
        ) const
        {
            return precondition(wT, rT, cmpt);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
