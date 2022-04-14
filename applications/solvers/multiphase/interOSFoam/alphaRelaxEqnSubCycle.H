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
    surfaceScalarField rhoPhiSum = 0.0*rhoPhicr;

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
#       include "alphaRelaxEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhicr;
    }

    rhoPhicr = rhoPhiSum;
}
else
{
#       include "alphaRelaxEqn.H"
}

interface.correct();

rho_cr == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
