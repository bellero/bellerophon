volScalarField rAU(1.0/UEqn().A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn().H();

if(runTime.outputTime())
{
    mesh.lookupObject<volVectorField>("U_0").write();

    UEqn().H()().write();
    volScalarField divNewPhi(fvc::div(mesh.Sf() & Uf));
    divNewPhi.rename("divNewPhi");
    divNewPhi.write();

    volScalarField divOldPhi(fvc::div(mesh.Sf() & Uf.oldTime()));
    divOldPhi.rename("divOldPhi");
    divOldPhi.write();

    surfaceScalarField phiHbyAUncorr
    (
        "phiHbyAUncorr",
        (fvc::interpolate(HbyA) & mesh.Sf())
    );
    volScalarField divPhiHbyAUncorr = fvc::div(phiHbyAUncorr);
    divPhiHbyAUncorr.rename("divPhiHbyAUncorr");
    divPhiHbyAUncorr.write();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, Uf)
);

if(runTime.outputTime())
{
    volScalarField divPhiHbyACorr = fvc::div(phiHbyA);
    divPhiHbyACorr.rename("divPhiHbyACorr1");
    divPhiHbyACorr.write();
    volScalarField divrAU("divrAU", fvc::div(fvc::interpolate(rAU)));
    divrAU.write();
    volScalarField divddtCorr("divddtCorr", fvc::div(fvc::ddtCorr(U, Uf)));
    divddtCorr.write();
    surfaceScalarField phiCorr("phiCorr", mesh.Sf() & (Uf.oldTime() - fvc::interpolate(U.oldTime())));
    phiCorr.write();
    surfaceScalarField phiOldTime("phiOldTime", mesh.Sf() & Uf.oldTime());
    phiOldTime.write();
    surfaceScalarField phiCorrUnlim("phiCorrUnlim", mag(phiCorr)/(mag(phiOldTime) + dimensionedScalar("small", phi.dimensions(), SMALL)));
    phiCorrUnlim.write();
    volScalarField divUOld("divUOld", fvc::div(mesh.Sf() & fvc::interpolate(U.oldTime())));
    divUOld.write();
    volScalarField divPhiCorr("divPhiCorr", fvc::div(mesh.Sf() & (Uf.oldTime() - fvc::interpolate(U.oldTime()))));
    divPhiCorr.write();

    volVectorField Uold("Uold", U.oldTime());
    Uold.write();

    surfaceVectorField Ufold("Ufold", Uf.oldTime());
    Ufold.write();
}

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent() && false)
{
    rAtU = 1.0/max(1.0/rAU - UEqn().H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

surfaceScalarField rAUf("rAUf", fvc::interpolate(rAtU()));

// Update the fixedFluxPressure BCs to ensure flux consistency
setSnGrad<fixedFluxPressureFvPatchScalarField>
(
    p.boundaryField(),
    (
        phiHbyA.boundaryField()
      - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())
    )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    volScalarField divHbyA = fvc::div(phiHbyA);

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == divHbyA
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
    if(runTime.outputTime())
        {
            volVectorField gradP("gradP", fvc::grad(p));
            gradP.write();
            HbyA.write();
            divHbyA.rename("divHbyA");
            divHbyA.write();
            rAtU().rename("rAtU");
            rAtU().write();
            volScalarField rauf = fvc::div(rAUf);
            rauf.rename("divrAUf");
            rauf.write();
            volScalarField pSol(p);
            pSol.rename("pSol");
            pSol.write();
            rAU.rename("rAU");
            rAU.write();
        }
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += n*(phi/mesh.magSf() - (n & Uf));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

if(runTime.outputTime())
    {
        volScalarField divUf("divUf", fvc::div(Uf & mesh.Sf()));
        divUf.write();
        volScalarField divUfold("divUfOld", fvc::div(Uf.oldTime() & mesh.Sf()));
        divUfold.write();
        volVectorField gradP("gradP_final", fvc::grad(p));
        gradP.write();
    }
