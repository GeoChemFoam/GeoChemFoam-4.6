{

    forAll(solutionSpecies, i)
	{
		volScalarField& Yi = speciesMixture.Y(i);
		dimensionedScalar DYi = speciesMixture.DY(i);

                dimensionedScalar Yimax = max(Yi);

		while (steadyState.correctNonOrthogonal())
		{
			fvScalarMatrix YiEqn
			(
				  fvm::div(phi, Yi, "div(phi,Yi)")
                                + fvm::SuSp(-fvc::div(phi),Yi)
				- fvm::laplacian(DYi, Yi)
			);

			YiEqn.relax();
			YiEqn.solve(mesh.solutionDict().solver("Yi"));
                        Yi=min(Yimax,max(Yi,0.0*Yi));
		}

		Info<< solutionSpecies[i] << " concentration = "
			<< "  Min(Yi) = " << gMin(Yi)
			<< "  Max(Yi) = " << gMax(Yi)
			<< endl;
	}	
}

speciesMixture.correct();
