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

#include "AxTx.H"
#include "bellerophonInterpolation.H"

#include "fvCFD.H"
#include "word.H"
#include "fvMesh.H"

#include "bellerophonPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AxTx, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<AxTx>
        addAxTxAsymMatrixConstructorToTable_;

    lduMatrix::solver::addsymMatrixConstructorToTable<AxTx>
        addAxTxSymMatrixConstructorToTable_;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AxTx::AxTx
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
            "Foam::AxTx::AxTx"
            "("
            "const word& fieldName, "
            "const lduMatrix& matrix, "
            "const FieldField<Field, scalar>& interfaceBouCoeffs, "
            "const FieldField<Field, scalar>& interfaceIntCoeffs, "
            "const lduInterfaceFieldPtrsList& interfaces, "
            "const dictionary& solverControls"
            ")"
        )<<"Using AxTx without bellerophonInterfaceField."<<endl;
    }

    bellerophon::Interpolation().updateMatrix(bellerophonMatrix_, firstInterface);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::AxTx::solve
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

    scalarField Ax(nCells);
    scalarField Tx(nCells);

    scalarField A(nCells*nCells);
    scalarField T(nCells*nCells);

    for(label cellI = 0; cellI < nCells; cellI++)
    {
        psi[cellI] = 1.0;
    }

    // --- Calculate A.psi2 and T.psi2
    bellerophonMatrix_.Amul(Ax, psi, interfaceBouCoeffs_, interfaces_, cmpt);
    bellerophonMatrix_.Tmul(Tx, psi, interfaceIntCoeffs_, interfaces_, cmpt);

    const fvMesh& mesh
    (
        bellerophonMatrix_.mesh().thisDb().lookupObject<fvMesh>("fvSchemes")
    );

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
    Ax_.internalField() = Ax;
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
    Tx_.internalField() = Tx;
    Tx_.write();


    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = 1;
    solverPerf.finalResidual() = 1;


    return solverPerf;
}


// ************************************************************************* //
