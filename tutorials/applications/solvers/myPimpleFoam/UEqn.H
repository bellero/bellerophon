// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

{
    volScalarField divPhiU = fvc::div(phi);
    divPhiU.rename("divPhiU");
    divPhiU.write();
}
tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

if(runTime.outputTime())
{
    volScalarField p0(p);
    p0.rename("p_U");
    p0.write();
    tmp<volVectorField> tGradP = fvc::grad(p);
    tGradP().rename("grapP_U");
    tGradP().write();
    volScalarField diag(p);
    diag.rename("diag");
    diag.internalField() = UEqn().diag();
    diag.write();
}

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
}
