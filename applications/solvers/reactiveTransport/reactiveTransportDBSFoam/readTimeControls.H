adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

maxDeltaEps =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaEps", 0.1);

maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

