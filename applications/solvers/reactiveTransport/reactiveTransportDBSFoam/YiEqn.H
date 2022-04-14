{    
    R = 0*R;
    forAll(kineticPhaseReactions,j)
    {
        volScalarField Rj = speciesMixture.Rate(j);

        forAll(kineticPhases,k)
        {
            if (kineticPhaseReactions[j]==kineticPhases[k])
            {
                volScalarField omega = speciesMixture.omega(k);
                Rj *= (1-omega);
            }
        }

        if (VoS=="VoS-psi")
        {
            R  += Mw/rhos*Rj*a;
        }
        else if (VoS=="iVoS")
        {
            surfaceScalarField phiR = Mw/rhos*fvc::interpolate(Rj)*nEpsf;

            surfaceScalarField phiREps =
                fvc::flux
                (
                    phiR,
                    eps,
                    "div(phiR,eps)"
                );


            R+=mag(eps*fvc::div(phiR)-fvc::div(phiREps))*sign(Rj);
        }
        else
        {
             Info << "VoS-scheme " << VoS << " does not exist" << endl;
             Info<< "VoS-scheme should be equal to VoS-psi or iVoS"
             << endl
             << abort(FatalError);
        } 
    }

    
    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = speciesMixture.Y(i);
        dimensionedScalar DYi = speciesMixture.DY(i);
	
        volScalarField Ri = 0*rhos/Mw*R;
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

            if (VoS=="VoS-psi")
            {
                Ri += scoeffj[i]*a*Rj; 
            }
            else if (VoS=="iVoS")
            {
                surfaceScalarField phiR = -scoeffj[i]*fvc::interpolate(Rj)*nEpsf;

                surfaceScalarField phiREps = fvc::flux(phiR,eps,"div(phiR,eps)");
                Ri += mag(fvc::div(phiREps)-eps*fvc::div(phiR))*sign(Rj);
            }
            else
            {
                Info<< "VoS-scheme should be equal to VoS-psi or iVoS"
                << endl
                << abort(FatalError);
            }
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
            Yi = max(Yi,0*Yi);
        }


        Info<< solutionSpecies[i] << " concentration = "
        << "  Min(Yi) = " << gMin(Yi)
        << "  Max(Yi) = " << gMax(Yi)
        << endl;
    }

    speciesMixture.correct();
}
