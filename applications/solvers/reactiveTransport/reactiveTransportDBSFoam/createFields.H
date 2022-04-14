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

    Info << "Reading field eps\n" << endl;
    volScalarField eps
    (
        IOobject
        (
            "eps",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    surfaceScalarField nEpsf
    (
       IOobject
       (
          "nEpsf",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
       ),
       mesh,
       dimensionedScalar("nEpsf",dimArea, 0.0)
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
    const wordList& kineticPhases = speciesMixture.kineticPhases();
    const wordList& kineticPhaseReactions = speciesMixture.kineticPhaseReactions();

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


	dimensionedScalar kf
	(
		transportProperties.lookup("kf")
	);

	dimensionedScalar Mw
	(
		transportProperties.lookup("Mw")
	);

    dimensionedScalar rhos
    (
        transportProperties.lookup("rhos")
    );

	//permeability
    volScalarField Kinv = kf*pow(1-eps,2)/pow(eps,3);

    Info << "Reading reaction rate if present\n" << endl;
    volScalarField R
    (
        IOobject
        (
            "R",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("R",dimless/dimTime, 0.0)
    );

    dimensionedScalar deltaN = 1e-13/pow(average(mesh.V()), 1.0/3.0);

    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i = 0; i < p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    word VoS(mesh.solutionDict().subDict("STEADYSTATE").lookup("VoS"));

    scalar psiCoeff(mesh.solutionDict().subDict("STEADYSTATE").lookupOrDefault("psiCoeff",4.0));

    bool adaptPsiCoeff(mesh.solutionDict().subDict("STEADYSTATE").lookupOrDefault("adaptPsiCoeff",false));

