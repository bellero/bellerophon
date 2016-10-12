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

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (simple.loop())
    {
        Info<<"Time = "<<runTime.timeName()<<nl<<endl;
        // --- Pressure-velocity SIMPLE corrector
        {
            //UEqn.H

            divPhi0=fvc::div(phi);

            // Momentum predictor
            tmp<fvVectorMatrix> UEqn
            (
                fvm::div(phi, U)
              - fvm::laplacian(turbulence->nuEff(), U)
              - fvc::div(turbulence->nuEff()*dev(T(fvc::grad(U))))
            //    ==
            //     sources(U)
            );

            deltaU=U;
            deltaU2=U;

            //UEqn().boundaryManipulate(U.boundaryField());

            gradP_U.internalField()=fvc::grad(p);

            rhs1.internalField() = fvm::div(phi, U)().source();

            UEqn().relax();

//             sources.constrain(UEqn());

            rhs2.internalField() = fvm::laplacian(turbulence->nuEff(), U)().source();

            rhs3.internalField() = fvc::div(turbulence->nuEff()*dev(T(fvc::grad(U))))().internalField();


            Info<<"Going to solve for U"<<endl;
            solve(UEqn()==-fvc::grad(p));

            deltaU-=U;

            U0=U;

            {
                p.boundaryField().updateCoeffs();

                UEqnDiag.internalField()=UEqn().diag();
                UEqnD.internalField()=UEqn().D();
                UEqnDbyV.internalField()=UEqnD/mesh.V();
                UEqnA=UEqn().A();
                rAU=1.0/UEqnA;

                UEqnH = UEqn().H()();
                UEqnLduH.internalField() = UEqn().lduMatrix::H(U.internalField());
                U = rAU*UEqnH;
                U1 = U;

                UEqn.clear();

//                 const scalarField oldPhi1 = phi.boundaryField()[3];
//                 const scalarField oldPhi2 = phi.boundaryField()[4];

                phi = fvc::interpolate(U, "interpolate(HbyA)") & mesh.Sf();

//                 phi.boundaryField()[3] = oldPhi1;
//                 phi.boundaryField()[4] = oldPhi2;

                phi1 = fvc::interpolate(U, "interpolate(HbyA)") & mesh.Sf();

                divPhi1=fvc::div(phi1);

//                adjustBellerophonPhi(phi, U, p);
                adjustPhi(phi, U, p);

                divPhi2=fvc::div(phi);

                phi2=phi;


                // Non-orthogonal pressure corrector loop
                while (simple.correctNonOrthogonal())
                {
                    divPhi=fvc::div(phi);

                    fvScalarMatrix laplacianRAuP(fvm::laplacian(rAU, p));

//                     runTime.write();

                    fvScalarMatrix pEqn
                    (
                        laplacianRAuP == divPhi
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    //pEqn.boundaryManipulate(p.boundaryField());

                    pEqn.solve();

                    if (simple.finalNonOrthogonalIter())
                    {
                        phi3 = pEqn.flux();
                        phi -= pEqn.flux();
                        divPhi3 = fvc::div(phi);
                        forAll(divPhi3.boundaryField(), patchI)
                        {
                            divPhi3.boundaryField()[patchI] = dynamic_cast<scalarField&>(phi.boundaryField()[patchI]);
                        }
                        gradP = fvc::grad(p);
                    }
                }

//                 #include "continuityErrs.H"
                {
                    contErr = fvc::div(phi);

                    scalar sumLocalContErr = runTime.deltaTValue()*
                        mag(contErr)().weightedAverage(mesh.V()).value();

                    scalar globalContErr = runTime.deltaTValue()*
                        contErr.weightedAverage(mesh.V()).value();
                    cumulativeContErr += globalContErr;

                    Info<< "time step continuity errors : sum local = " << sumLocalContErr
                        << ", global = " << globalContErr
                        << ", cumulative = " << cumulativeContErr
                        << endl;

                }

                // Explicitly relax pressure for momentum corrector
                unrelaxedP = p;
                p.relax();


                // Momentum corrector
                U -= rAU*fvc::grad(p);

                gradP = fvc::grad(p);
            }
        }

        deltaU2-=U;

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
