{
	//Transport step
	Info << "transport step" << endl;

    forAll(solutionSpecies, i)
	{
		volScalarField& Yi = speciesMixture.Y(i);
		dimensionedScalar DYi = speciesMixture.DY(i);

		while (simple.correctNonOrthogonal())
		{
			fvScalarMatrix YiEqn
			(
				fvm::ddt(Yi)
				+ fvm::div(phi, Yi, "div(phi,Yi)")
				- fvm::laplacian(DYi, Yi)
			);

			//solve equations
		    if (simple.finalNonOrthogonalIter())
		    {
		        YiEqn.solve(mesh.solutionDict().solver("YiFinal"));
		    }
		    else
		    {
				YiEqn.solve(mesh.solutionDict().solver("Yi"));
			}
		}

		Info<< solutionSpecies[i] << " concentration = "
			<< "  Min(Yi) = " << gMin(Yi)
			<< "  Max(Yi) = " << gMax(Yi)
			<< endl;
	}

    //Correct reaction thermo
    Info << "reaction thermo step" << endl;
    speciesMixture.correct();
}
