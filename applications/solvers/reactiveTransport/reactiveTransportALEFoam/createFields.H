    Info << "Reading field p\n" << endl;
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

    Info << "Reading field U\n" << endl;
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

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

#   include "createPhi.H"

	Info << "Create species mixture\n" << endl;
    multiComponentTransportMixture<reactiveMixture> speciesMixture(mesh);
    const speciesTable& solutionSpecies = speciesMixture.species();

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());

    singlePhaseTransportModel laminarTransport(U, phi);

    dimensionedScalar nu
    (
            transportProperties.lookup("nu")
    );

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

	label checkFrequency =  	runTime.controlDict().lookupOrDefault<label>("checkFrequency", 1.0);

        label remeshFrequency =         runTime.controlDict().lookupOrDefault<label>("remeshFrequency", 10000);


