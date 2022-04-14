    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field pvd\n" << endl;
    volScalarField pvd
    (
        IOobject
        (
            "pvd",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pd
    );

    Info<< "Reading field pcr\n" << endl;
    volScalarField pcr
    (
        IOobject
        (
            "pcr",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pd
    );

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
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

	//viscous drag velocity
    volVectorField Uvd
    (
        IOobject
        (
            "Uvd",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U
    );

	//capillary relaxation velocity
    Info<< "Reading field Ucr\n" << endl;
    volVectorField Ucr
    (
        IOobject
        (
            "Ucr",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
		U
    );

	Ucr.boundaryField() == 0*U.boundaryField();

#   include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "alpha1");

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    // Need to store rho_vd for ddt(rho, U_vd)
    volScalarField rho_vd
    (
        IOobject
        (
            "rho_vd",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho_vd.oldTime();

    // Need to store rho_cr for ddt(rho, U_cr)
    volScalarField rho_cr
    (
        IOobject
        (
            "rho_cr",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho_cr.oldTime();

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );

    surfaceScalarField rhoPhivd
    (
        IOobject
        (
            "rho*phivd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phivd
    );

    surfaceScalarField rhoPhicr
    (
        IOobject
        (
            "rho*phicr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phicr
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("gh", g & mesh.Cf());

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

	const volScalarField& pc = interface.pc();
    mesh.schemesDict().setFluxRequired(pc.name());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd + pc + rho_vd*gh
    );


    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, pimple.dict(), pdRefCell, pdRefValue);
    mesh.schemesDict().setFluxRequired(pvd.name());
    mesh.schemesDict().setFluxRequired(pcr.name());

    scalar pRefValue = 0.0;

    if (pd.needReference())
    {
        pRefValue = readScalar(pimple.dict().lookup("pRefValue"));

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pdRefCell)
        );
    }

     // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    scalar deltaTCR = readScalar(pimple.dict().lookup("deltaTCR"));
    label relaxMin = readLabel(pimple.dict().lookup("relaxMin"));

    scalar residual = readScalar(pimple.dict().lookup("residual"));

    scalar nRelaxSteps =0;
    scalar nSteps =0;
    scalar nRelaxTot = 0;
    label nIterLast = 1;
    
