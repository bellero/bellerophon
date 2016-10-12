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
    pimpleFoam

Description
    Large time-step transient solver for incompressible, flow using the PIMPLE
    (merged PISO-SIMPLE) algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable finite volume options, e.g. MRF, explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // Solve the Momentum equation

                MRF.correctBoundaryVelocity(U);

                divPhiU = fvc::div(phi);

                tmp<fvVectorMatrix> UEqn
                (
                    fvm::ddt(U)
                    + fvm::div(phi, U)
                    + MRF.DDt(U)
                    + turbulence->divDevReff(U)
                   ==
                    fvOptions(U)
                );

                UEqn().relax();

                gradP = fvc::grad(p);

                fvOptions.constrain(UEqn());

                UEqnDiag.internalField()=UEqn().D();

                if (pimple.momentumPredictor())
                {
                    solve(UEqn() == -fvc::grad(p));

                    fvOptions.correct(U);
                }

                U0=U;


            // --- Pressure corrector loop
            while (pimple.correct())
            {
                // Solve pressure correction
                {
                    rAU = 1.0/UEqn().A();

                    HbyA = rAU*UEqn().H();

                    phiHbyA = (fvc::interpolate(HbyA) & mesh.Sf())
                            + fvc::interpolate(rAU)*fvc::ddtCorr(U,phi);

                    MRF.makeRelative(phiHbyA);

                    adjustPhi(phiHbyA, U, p);

                    rAtU=rAU;

                    if(pimple.consistent())
                    {
                        rAtU = 1.0/max(1.0/rAU - UEqn().H1(), 0.1/rAU);
                        phiHbyAcorr = fvc::interpolate(rAtU - rAU)*fvc::snGrad(p)*mesh.magSf();
                        phiHbyA += phiHbyAcorr;
                        HbyAcorr = (rAU - rAtU)*gradP;
                        HbyA -= HbyAcorr;
                    }

                    if (pimple.nCorrPISO() <= 1)
                    {
                        UEqn.clear();
                    }

                    rAUf = fvc::interpolate(rAU);

                    // Update the fixedFluxPressure BCs to ensure flux consistency
                    setSnGrad<fixedFluxPressureFvPatchScalarField>
                    (
                        p.boundaryField(),
                        (
                            phiHbyA.boundaryField()
                          - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())
                        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
                    );

                    while (pimple.correctNonOrthogonal())
                    {
                        divPhi = fvc::div(phiHbyA);

                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAUf, p) == divPhi
                        );

                        pEqn.setReference(pRefCell, pRefValue);

                        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

                        if (pimple.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }

                    #include "continuityErrs.H"

                    // Explicitly relax pressure for momentum corrector
                    p.relax();

                    gradPCorr = fvc::grad(p);
                    Ucorr = rAtU * fvc::grad(p);
                    U = HbyA - Ucorr;

                    U.correctBoundaryConditions();
                    fvOptions.correct(U);
                }

            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
