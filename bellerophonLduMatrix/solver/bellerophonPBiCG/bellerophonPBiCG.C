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

#include "bellerophonPBiCG.H"
#include "bellerophonInterpolation.H"

#include "fvCFD.H"
#include "word.H"
#include "fvMesh.H"

#include "bellerophonPreconditioner.H"
#include "printMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonPBiCG, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<bellerophonPBiCG>
        addbellerophonPBiCGAsymMatrixConstructorToTable_;

    lduMatrix::solver::addsymMatrixConstructorToTable<bellerophonPBiCG>
        addbellerophonPBiCGSymMatrixConstructorToTable_;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bellerophonPBiCG::bellerophonPBiCG
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
            "Foam::bellerophonPBiCG::bellerophonPBiCG"
            "("
            "const word& fieldName, "
            "const lduMatrix& matrix, "
            "const FieldField<Field, scalar>& interfaceBouCoeffs, "
            "const FieldField<Field, scalar>& interfaceIntCoeffs, "
            "const lduInterfaceFieldPtrsList& interfaces, "
            "const dictionary& solverControls"
            ")"
        )<<"Using bellerophonPBiCG without bellerophonInterfaceField."<<endl;
    }

    const_cast<bellerophonInterfaceField*>(firstInterface)->initEvaluate();

    bellerophon::Interpolation().updateMatrix(bellerophonMatrix_, firstInterface);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::bellerophonPBiCG::solve
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

    register label nCells = psi.size();

    scalar* __restrict__ psiPtr = psi.begin();

    scalarField pA(nCells);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField pT(nCells, 0.0);
    scalar* __restrict__ pTPtr = pT.begin();

    scalarField wA(nCells);
    scalar* __restrict__ wAPtr = wA.begin();

    scalarField wT(nCells);
    scalar* __restrict__ wTPtr = wT.begin();

    scalar wArT = solverPerf.great_;
    scalar wArTold = wArT;

    // --- Calculate A.psi and T.psi
    bellerophonMatrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);
    bellerophonMatrix_.Tmul(wT, psi, interfaceIntCoeffs_, interfaces_, cmpt);

    tmp<scalarField> tCorrSource =
        bellerophon::Interpolation().correctSource
        (
            source, bellerophonMatrix_.diag(), cmpt
        );

    const scalarField& source2 = tCorrSource();

    // --- Calculate initial residual and transpose residual fields
    scalarField rA(source2 - wA);
    scalarField rT(source2 - wT);

    if(debug)
    {
        const fvMesh& mesh
        (
            bellerophonMatrix_.mesh().thisDb().lookupObject<fvMesh>("fvSchemes")
        );

        volScalarField rA_
        (
            IOobject
            (
                "rA_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        rA_.internalField() = rA;
        rA_.write();

        volScalarField rT_
        (
            IOobject
            (
                "rT_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        rT_.internalField() = rT;
        rT_.write();

        volScalarField source_
        (
            IOobject
            (
                "source_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        source_.internalField() = source2;
        source_.write();

        volScalarField psi_
        (
            IOobject
            (
                "psi_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        psi_.internalField() = psi;
        psi_.write();

    }

    scalar* __restrict__ rAPtr = rA.begin();
    scalar* __restrict__ rTPtr = rT.begin();

    // --- Calculate normalisation factor
    scalar normFactor = this->normFactor(psi, source2, wA, pA);

    if (bellerophonLduMatrix::debug >= 2 || debug)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if (!solverPerf.checkConvergence(tolerance_, relTol_))
    {
        // --- Select and construct the preconditioner
        autoPtr<bellerophonLduMatrix::preconditioner> preconPtr =
        bellerophonLduMatrix::preconditioner::New
        (
            bellerophonMatrix_,
            controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Store previous wArT
            wArTold = wArT;

            // --- Precondition residuals
            preconPtr->precondition(wA, rA, cmpt);
            preconPtr->preconditionT(wT, rT, cmpt);

            // --- Update search directions:
            wArT = gSumProd(wA, rT);

            if (solverPerf.nIterations() == 0)
            {
                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                    pTPtr[cell] = wTPtr[cell];
                }
            }
            else
            {
                if(wArTold == 0.0) wArTold = SMALL;
                scalar beta = wArT/wArTold;

                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                    pTPtr[cell] = wTPtr[cell] + beta*pTPtr[cell];
                }
            }


            // --- Update preconditioned residuals
            bellerophonMatrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);
            bellerophonMatrix_.Tmul(wT, pT, interfaceIntCoeffs_, interfaces_, cmpt);

            scalar wApT = gSumProd(wA, pT);


            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApT)/normFactor)) break;

            // --- Update solution and residual:

            scalar alpha = wArT/wApT;

            for (register label cell=0; cell<nCells; cell++)
            {
                psiPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
                rTPtr[cell] -= alpha*wTPtr[cell];
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;

        } while
        (
            (
                solverPerf.nIterations()++ < maxIter_
             && !(solverPerf.checkConvergence(tolerance_, relTol_))
            )
            ||
            solverPerf.nIterations() < minIter_
        );
    }

    //DEBUGGING OUTPUT
    if(debug)
    {
        const fvMesh& mesh
        (
            bellerophonMatrix_.mesh().thisDb().lookupObject<fvMesh>("fvSchemes")
        );

        volScalarField result
        (
            IOobject
            (
                "result_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        result.internalField() = psi;
        result.write();

        volScalarField pA_final
        (
            IOobject
            (
                "pA_final_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        pA_final.internalField() = pA;
        pA_final.write();

        volScalarField pT_final
        (
            IOobject
            (
                "pT_final_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        pT_final.internalField() = pT;
        pT_final.write();

        volScalarField wA_final
        (
            IOobject
            (
                "wA_final_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        wA_final.internalField() = wA;
        wA_final.write();

        volScalarField wT_final
        (
            IOobject
            (
                "wT_final_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        wT_final.internalField() = wT;
        wT_final.write();

        volScalarField rA_final
        (
            IOobject
            (
                "rA_final_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        rA_final.internalField() = rA;
        rA_final.write();

        volScalarField rT_final
        (
            IOobject
            (
                "rT_final_"+fieldName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
        );
        rT_final.internalField() = rT;
        rT_final.write();

//         volScalarField sourceField
//         (
//             IOobject
//             (
//                 "source",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("zero",dimless,0.0)
//         );
//         sourceField.internalField() = source2;
//         sourceField.write();
    }

    return solverPerf;
}


// ************************************************************************* //
