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

\*---------------------------------------------------------------------------*/

#include "printMatrix.H"
#include "bellerophonInterpolation.H"

#include "fvCFD.H"
#include "word.H"
#include "fvMesh.H"

#include "bellerophonPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(printMatrix, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<printMatrix>
        addprintMatrixAsymMatrixConstructorToTable_;

    lduMatrix::solver::addsymMatrixConstructorToTable<printMatrix>
        addprintMatrixSymMatrixConstructorToTable_;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::printMatrix::printMatrix
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    bellerophonMatrix_(matrix)
{
    label firstBellerophonInterface = -1;
    forAll(interfaces,interfaceI)
    {
        if
        (
            firstBellerophonInterface<0
            &&
            interfaces.set(interfaceI)
            &&
            isA<bellerophonInterfaceField>(interfaces[interfaceI])
        )
        {
            firstBellerophonInterface=interfaceI;
        }
    }

    bellerophonInterfaceField const * firstInterface(NULL);

    if(firstBellerophonInterface >= 0)
    {
        firstInterface = dynamic_cast< bellerophonInterfaceField const * >
            (interfaces(firstBellerophonInterface));
    }
    else
    {
        WarningIn
        (
            "Foam::printMatrix::printMatrix"
            "("
            "const word& fieldName, "
            "const lduMatrix& matrix, "
            "const FieldField<Field, scalar>& interfaceBouCoeffs, "
            "const FieldField<Field, scalar>& interfaceIntCoeffs, "
            "const lduInterfaceFieldPtrsList& interfaces, "
            "const dictionary& solverControls"
            ")"
        )<<"Using printMatrix without bellerophonInterfaceField."<<endl;
    }

    bellerophon::Interpolation().updateMatrix(bellerophonMatrix_, firstInterface);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::printMatrix::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        bellerophonLduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

//     bellerophonMatrix_.prepareMultiplications(psi);

    register label nCells = psi.size();
    scalarField psi2(psi);

    scalarField Ax(nCells);
    scalarField Tx(nCells);

    scalarField A(nCells*nCells);
    scalarField T(nCells*nCells);

    for(label cellI = 0; cellI < nCells; cellI++)
    {
        forAll(psi2, psi2I)
        {
            if(psi2I == cellI)
            {
                psi2[psi2I] = 1.0;
            }
            else
            {
                psi2[psi2I] = 0.0;
            }
        }
        // --- Calculate A.psi2 and T.psi2
        bellerophonMatrix_.Amul(Ax, psi2, interfaceBouCoeffs_, interfaces_, cmpt);
        bellerophonMatrix_.Tmul(Tx, psi2, interfaceIntCoeffs_, interfaces_, cmpt);

//         Pout<<"psi2:"<<psi2<<nl<<"Ax: "<<Ax<<nl<<"Tx: "<<Tx<<endl;

        forAll(Ax, AxI)
        {
            A[AxI*nCells+cellI]=Ax[AxI];
            T[AxI*nCells+cellI]=Tx[AxI];
        }
    }

    const fvMesh& mesh
    (
        bellerophonMatrix_.mesh().thisDb().lookupObject<fvMesh>("fvSchemes")
    );

    word prefix;
    if(Pstream::parRun())
    {
        prefix=mesh.time().timeName()+"_"+name(Pstream::myProcNo())+"_";
    }
    else
    {
        prefix=mesh.time().timeName()+"_";
    }


    OFstream Afile(prefix+"A.dat");
    OFstream Tfile(prefix+"T.dat");

    for(label cellI = 0; cellI < nCells; cellI++)
    {
        const label row = cellI*nCells;
        Afile<<A[row];
        Tfile<<T[row];
        for(label cellII = 1; cellII < nCells; cellII++)
        {
            Afile<<"\t"<<A[row+cellII];
            Tfile<<"\t"<<T[row+cellII];
        }
        Afile<<endl;
        Tfile<<endl;
    }

    if(debug)
    {

        volScalarField Ax_
        (
            IOobject
            (
                "Ax_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        Ax_.primitiveFieldRef() = Ax;
        Ax_.write();

        volScalarField Tx_
        (
            IOobject
            (
                "Tx_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        Tx_.primitiveFieldRef() = Tx;
        Tx_.write();

    }


    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = 1;
    solverPerf.finalResidual() = 1;


    return solverPerf;
}


// ************************************************************************* //
