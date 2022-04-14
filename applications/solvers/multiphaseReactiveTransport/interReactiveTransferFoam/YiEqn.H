{
    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = speciesMixture.Y(i);

        //compute diffusion coefficient
        surfaceScalarField D = speciesMixture.DmY(i);

        //calculate Henry's transfer flux
        surfaceScalarField phiH = speciesMixture.phiH(i);
        surfaceScalarField phiHUp = speciesMixture.phiHUp(i);
        surfaceScalarField phiHDown = speciesMixture.phiHDown(i);

        volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

        while (simple.correctNonOrthogonal())
        {
            //compute Y equation
            fvScalarMatrix YiEqn
            (
                fvm::ddt(Yi)
                - fvm::laplacian(D,Yi)
                ==
                SpY
            );

            if (phiHScheme == "Gauss linear") YiEqn += fvm::div(phiH,Yi,"div(phiH,Yi)");
            else if (phiHScheme == "Gauss upwind")
            {
                YiEqn += fvm::div(phiHUp,Yi,"div(phiH,Yi)")
                + fvm::div(phiHDown,Yi,"div(phiH,Yi)");
            }
            else
            {
                Info<< "div(phiH,Yi) should be equal to Gauss linear or Gauss upwind"
				<< endl
				<< abort(FatalError);
            }

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

        Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
        << "  Min(Yi) = " << gMin(Yi.internalField())
        << "  Max(Yi) = " << gMax(Yi.internalField())
        << endl;
    }

    //calculate mass transfer flux
    speciesMixture.correct();
 }
