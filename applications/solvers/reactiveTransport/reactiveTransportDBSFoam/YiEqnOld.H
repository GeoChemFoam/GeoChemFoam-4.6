{    

    volScalarField a = mag(fvc::grad(eps))*(4*eps*(1-eps));

    volVectorField gradEps = fvc::grad(eps);
    surfaceVectorField gradEpsf = fvc::interpolate(gradEps);
    surfaceVectorField nEpsv = gradEpsf/(mag(gradEpsf) + deltaN);
    surfaceScalarField nEpsf = nEpsv & mesh.Sf();


    R = 0*R;
    
    forAll(solutionSpecies, i)
	{
		volScalarField& Yi = speciesMixture.Y(i);
		dimensionedScalar DYi = speciesMixture.DY(i);
	
                volScalarField Ri = 0*R*rhos/Mw;
		
		forAll(kineticPhaseReactions,j)
		{
		      volScalarField Rj = speciesMixture.Rate(j);
                      const scalarList& scoeffj = speciesMixture.scoeff(j);
        
		       forAll(kineticPhases,k)
		       {
		           if (kineticPhaseReactions[j]==kineticPhases[k])
                          {
		               volScalarField omega = speciesMixture.omega(k);
                               Rj *= (1-omega);
		           }
		       }

                       Ri += scoeffj[i]*Rj*a;
                       R  += Mw/rhos*Rj*a;
		}
    
                volScalarField D = eps*DYi;
        
		while (steadyState.correctNonOrthogonal())
		{
			fvScalarMatrix YiEqn
			(
				  fvm::div(phi, Yi, "div(phi,Yi)")
                                + fvm::SuSp(-fvc::div(phi),Yi)
				- fvm::laplacian(D, Yi)
				+ Ri 
			);

			YiEqn.relax();
			YiEqn.solve(mesh.solutionDict().solver("Yi"));
		}


		Info<< solutionSpecies[i] << " concentration = "
			<< "  Min(Yi) = " << gMin(Yi)
			<< "  Max(Yi) = " << gMax(Yi)
			<< endl;
	}

	speciesMixture.correct();
}
