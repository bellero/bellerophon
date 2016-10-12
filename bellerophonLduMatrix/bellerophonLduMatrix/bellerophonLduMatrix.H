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
    Foam::bellerophonLduMatrix

Description
    bellerophonLduMatrix is a extended version of the standard lduMatrix - you
    might have guessed it. It extends the matrix with off-diagnoal coefficients
    for the implicit coupling on overset grids as well as a modified copy of
    the off-diagonal coefficients.

    Objects of this class are initialised after selection of the bellerophon
    solver from the fvSolution dictionary.
    The solver instances a matrix of this class and
    (e.g. bellerophonPBiCG), so no new solvers are derived from this class.

\*---------------------------------------------------------------------------*/

#ifndef bellerophonLduMatrix_H
#define bellerophonLduMatrix_H

#include "lduMatrix.H"
#include "lduMesh.H"

#include "interpolationItem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class bellerophonLduMatrix;
Ostream& operator<<(Ostream&, const bellerophonLduMatrix&);


/*---------------------------------------------------------------------------*\
                           Class bellerophonLduMatrix Declaration
\*---------------------------------------------------------------------------*/

class bellerophonLduMatrix
:
public lduMatrix
{
    // private data

        //- Modified matrix coefficients
        scalarField *modUpperPtr_, *modLowerPtr_;

        //- List of primary Donor cells
        const labelList& primaryDonorCells_;

        //- Interpolation weights for primary donor cells
        const scalarField& primaryDonorWeights_;

        //- Interpolation items on own processor
        const List<interpolationItem>& ownInterpolationItems_;

        //- Interpolation items for other processors
        const List< List<interpolationItem> >& neighbourInterpolationItems_;

        const labelListList& neighbourValueToFieldMap_;

        mutable scalarField interpolatedPsi_;

        //- Buffer for acceptor values
        mutable scalarListList acceptorRowBuf_;

        //- Buffer for donor values
        mutable scalarListList donorColBuf_;

        //- Index of the donor on the local domain stored separatly for
        //  each domain the acceptor cell is in
        //    Amul: Holds the index of interpolated Psi which needs to be send
        //          to the corresponding proc
        //          (transfer raw value, needs to be multiplied by diagnoal
        //          element for that row)
        //    Tmul: Holds the index of the row the received value from the
        //          corresponding proc needs to be added to in interpolated Psi
        //          (transfer value mutilplied with diagonal coeff,
        //          directly added to that row )
        const labelListList& donorCols_;

        //- Rows of the acceptor cells on the local domain stored separatly for
        //  each domain the donor cell is in
        //    Amul: Holds the index of the row the received value from the
        //          corresponding proc needs to be added to
        //          (raw value transfered, needs to be multiplied by diagnoal
        //          element for that row)
        //    Tmul: Holds the index of Psi which needs to be send to the
        //          corresponding proc
        //          (transfer value mutilplied with diagonal coeff,
        //          directly added to that row )
        const labelListList& acceptorRows_;

        //- Buffer for interpolated psi values for neighbour processors
        mutable scalarListList neighbourInterpolationBuf_;

        //- Buffer for interpolated psi values from neighbour processors
        mutable scalarListList neighbourValueBuf_;

        //- Resize buffer if required
        void resizeBufs() const;

public:

/*    //- Abstract base-class for lduMatrix solvers
    //  obsolete, since apropriate solvers (e.g. bellerophonPBiCG) are selected
    //  before bellerophonLduMatrices are created
    class solver
    {
    protected:

        // Protected data

            word fieldName_;
            const lduMatrix& matrix_;
            const FieldField<Field, scalar>& interfaceBouCoeffs_;
            const FieldField<Field, scalar>& interfaceIntCoeffs_;
            lduInterfaceFieldPtrsList interfaces_;

            //- dictionary of controls
            dictionary controlDict_;

            //- Maximum number of iterations in the solver
            label maxIter_;

            //- Final convergence tolerance
            scalar tolerance_;

            //- Convergence tolerance relative to the initial
            scalar relTol_;


        // Protected Member Functions

            //- Read the control parameters from the controlDict_
            virtual void readControls();


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                symMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& interfaceBouCoeffs,
                    const FieldField<Field, scalar>& interfaceIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& solverControls
                ),
                (
                    fieldName,
                    matrix,
                    interfaceBouCoeffs,
                    interfaceIntCoeffs,
                    interfaces,
                    solverControls
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                bellerophonMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const FieldField<Field, scalar>& interfaceBouCoeffs,
                    const FieldField<Field, scalar>& interfaceIntCoeffs,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& solverControls
                ),
                (
                    fieldName,
                    matrix,
                    interfaceBouCoeffs,
                    interfaceIntCoeffs,
                    interfaces,
                    solverControls
                )
            );


        // Constructors

            solver
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& interfaceBouCoeffs,
                const FieldField<Field, scalar>& interfaceIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& solverControls
            );

        // Selectors

            //- Return a new solver
            static autoPtr<solver> New
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& interfaceBouCoeffs,
                const FieldField<Field, scalar>& interfaceIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& solverControls
            );



        //- Destructor
        virtual ~solver()
        {}


        // Member functions

            // Access

                const word& fieldName() const
                {
                    return fieldName_;
                }

                const lduMatrix& matrix() const
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


            //- Read and reset the solver parameters from the given stream
            virtual void read(const dictionary&);

            virtual solverPerformance solve
            (
                scalarField& psi,
                const scalarField& source,
                const direction cmpt=0
            ) const = 0;

            //- Return the matrix norm used to normalise the residual for the
            //  stopping criterion
            scalar normFactor
            (
                const scalarField& psi,
                const scalarField& source,
                const scalarField& Apsi,
                scalarField& tmpField
            ) const;
    };
*/


    class smoother;



   //- Abstract base-class for bellerophonLduMatrices preconditioners
    class preconditioner;

    // Static data

        // Declare name of the class and its debug switch
        ClassName("bellerophonLduMatrix");

        //- Large scalar for the use in solvers
        static const scalar great_;

        //- Small scalar for the use in solvers
        static const scalar small_;

    // Constructors

        //- Construct as copy
        bellerophonLduMatrix(const bellerophonLduMatrix&);


        //- Construct from lduMatrix
        bellerophonLduMatrix(const lduMatrix& );


    //- Destructor
        ~bellerophonLduMatrix();


    // Member functions

        // Access to coefficients

            //- Constant access to the modified lower triangle of the matrix
            inline const scalarField& modLower() const;

            //- Non-constant access to the modified lower triangle of the matrix
            inline scalarField& modLower();

            //- Constant access to the modified upper triangle of the matrix
            inline const scalarField& modUpper() const;

            //- Non-constant access to the modified upper triangle of the matrix
            inline scalarField& modUpper();

            //- Local phi to proc addressing for Amul
            inline const labelListList& donorCol() const;

            //- Proc to local Phi addressing for Amul
            inline const labelListList& acceptorRow() const;

        // Operations

            //- Sum off-diag elements to diag
            void sumDiag();

            //- Negative sum off-diag elements to diag
            void negSumDiag();

            //- Sum of magnitudes of off-diag elements
            void sumMagOffDiag(scalarField& sumOff) const;

            //- Matrix multiplication with updated interfaces.
            void Amul
            (
                scalarField&,
                const tmp<scalarField>&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;

            //- Matrix transpose multiplication with updated interfaces.
            void Tmul
            (
                scalarField&,
                const tmp<scalarField>&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;


            //- Sum the coefficients on each row of the matrix
            void sumA
            (
                scalarField&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&
            ) const;


            void residual
            (
                scalarField& rA,
                const scalarField& psi,
                const scalarField& source,
                const FieldField<Field, scalar>& interfaceBouCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;

            tmp<scalarField> residual
            (
                const scalarField& psi,
                const scalarField& source,
                const FieldField<Field, scalar>& interfaceBouCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;


            template<class Type>
            tmp<Field<Type> > H(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type> > H(const tmp<Field<Type> >&) const;

            tmp<scalarField> H1() const;

            template<class Type>
            tmp<Field<Type> > faceH(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type> > faceH(const tmp<Field<Type> >&) const;


    // Member operators

        void operator=(const bellerophonLduMatrix&);

        void negate();

        void operator+=(const bellerophonLduMatrix&);
        void operator-=(const bellerophonLduMatrix&);

        void operator*=(const scalarField&);
        void operator*=(scalar);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const bellerophonLduMatrix&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "bellerophonLduMatrixI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
