label nAlphaCorr
(
    readLabel(pimple.dict().lookup("nAlphaCorr"))
);

label nAlphaSubCycles
(
    readLabel(pimple.dict().lookup("nAlphaSubCycles"))
);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum = 0.0*rhoPhivd;

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
#       include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhivd;
    }

    rhoPhivd = rhoPhiSum;
}
else
{
#       include "alphaEqn.H"
}

//interface.correct();

rho_vd == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
