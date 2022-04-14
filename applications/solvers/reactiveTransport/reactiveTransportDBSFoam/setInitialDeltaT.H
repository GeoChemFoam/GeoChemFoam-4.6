
if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && (deltaEpsMax > maxDeltaEps))
    {
        Info << "Hello" << endl;
        runTime.setDeltaT
        (
            min
            (
                maxDeltaEps*runTime.deltaT().value()/deltaEpsMax,
                maxDeltaT
            )
        );
    }
}

// ************************************************************************* //
