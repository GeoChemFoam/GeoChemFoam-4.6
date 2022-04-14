
if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && (meshCoNum > maxMeshCo))
    {
        runTime.setDeltaT
        (
            min
            (
                maxMeshCo*runTime.deltaT().value()/meshCoNum,
                maxDeltaT
            )
        );
    }
}

// ************************************************************************* //
