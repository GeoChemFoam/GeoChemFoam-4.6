{
    word alpharScheme("div(phirb,alpha)");
    word YiScheme("div(phi,Yi)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(cYi*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = speciesMixture.Y(i);

        scalar maxYi = max(gMax(Yi),gMax(Yi.boundaryField()));

        Yi.oldTime() == Yi.oldTime()/maxYi;
        Yi == Yi/maxYi;

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

        //interface compression flux
        surfaceScalarField phiY =
        fvc::flux
        (
            phi,
            Yi,
            YiScheme
        )
        +phiComp*speciesMixture.compressionCoeff(i);

        MULES::explicitSolve(Yi, phi, phiY, 1.0, 0);

        Yi.oldTime() == Yi.oldTime()*maxYi;
        Yi == Yi*maxYi;
    }
 }
