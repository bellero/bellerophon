    Info<< "Reading field p" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField unrelaxedP
    (
        IOobject
        (
            "unrelaxedP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p
    );

    Info<< "Reading field U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U0
    (
        IOobject
        (
            "U0",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
     U
    );

    volVectorField U1
    (
        IOobject
        (
            "U1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
     U
    );

    volTensorField gradU
    (
        IOobject
        (
            "gradU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(U)
    );

    volVectorField deltaU
    (
        IOobject
        (
            "deltaU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
     U
    );

    volVectorField deltaU2
    (
        IOobject
        (
            "deltaU2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
     U
    );

    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,1,0,0,0,0), 0.0)
    );

    volVectorField UEqnH
    (
        IOobject
        (
            "UEqnH",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
     mesh,
     dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    );

    volVectorField UEqnLduH
    (
        IOobject
        (
            "UEqnLdiH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    );

    // main diagonal
    volScalarField UEqnDiag
    (
        IOobject
        (
            "UEqnDiag",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volVectorField laplacianNuU
    (
        IOobject
        (
            "laplacianNuU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    );

    volVectorField divNuEffDevTU
    (
        IOobject
        (
            "divNuEffDevTU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    );

    // main diagonal
    volScalarField UEqnD
    (
        IOobject
        (
            "UEqnD",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volScalarField UEqnDbyV
    (
        IOobject
        (
            "UEqnDbyV",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volScalarField V
    (
        IOobject
        (
            "V",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0), 0.0)
    );
    V.internalField()=mesh.V();

    // coefficients on the main diagonal
    volScalarField UEqnA
    (
        IOobject
        (
            "UEqnA",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    volScalarField contErr
    (
        IOobject
        (
            "contErr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    // div(phi)
    volScalarField divPhi
    (
        IOobject
        (
            "divPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    volScalarField divPhi0
    (
        IOobject
        (
            "divPhi0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    volScalarField divPhi1
    (
        IOobject
        (
            "divPhi1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    volScalarField divPhi2
    (
        IOobject
        (
            "divPhi2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    volScalarField divPhi3
    (
        IOobject
        (
            "divPhi3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    volScalarField divPhi4
    (
        IOobject
        (
            "divPhi4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    volVectorField gradP_U
    (
        IOobject
        (
            "gradP_U",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    );

    volVectorField gradP
    (
        IOobject
        (
            "gradP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    );

    volScalarField deltaP
    (
        IOobject
        (
            "deltaP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p
    );

    volVectorField rhs1
    (
        IOobject
        (
            "rhs1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimless, vector::zero)
    );

    volVectorField rhs2
    (
        IOobject
        (
            "rhs2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimless, vector::zero)
    );

    volVectorField rhs3
    (
        IOobject
        (
            "rhs3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimless, vector::zero)
    );


    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,3,-1,0,0,0,0), 0.0)
    );

    // second phi, from solved U field
    surfaceScalarField phi1
    (
        IOobject
        (
            "phi1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,3,-1,0,0,0,0), 0.0)
    );

    // third phi, after adjusting
    surfaceScalarField phi2
    (
        IOobject
        (
            "phi2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,3,-1,0,0,0,0), 0.0)
    );

    // fourth phi, from final U field
    surfaceScalarField phi3
    (
        IOobject
        (
            "phi3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,3,-1,0,0,0,0), 0.0)
    );

    // to faces interpolated U
    surfaceScalarField phi4
    (
        IOobject
        (
            "phi4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,3,-1,0,0,0,0), 0.0)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);

    Info<<"Creating turbulence"<<endl;

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
