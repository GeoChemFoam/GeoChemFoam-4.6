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


    multiComponentTransportMixture<phreeqcMixture> speciesMixture(mesh);
    const speciesTable& solutionSpecies = speciesMixture.species();
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	

