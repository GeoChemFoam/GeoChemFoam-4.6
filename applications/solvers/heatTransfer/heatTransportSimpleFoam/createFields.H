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

#   include "createPhi.H"


    Info<< "Reading field T\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField eps
    (
        IOobject
        (
           "eps",
           runTime.timeName(),
           mesh,
           IOobject::READ_IF_PRESENT,
           IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("eps",dimless,1.0)
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


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DTf
    (
        transportProperties.lookup("DTf")
    );

    dimensionedScalar DTs
    (
        transportProperties.lookup("DTs")
    );
