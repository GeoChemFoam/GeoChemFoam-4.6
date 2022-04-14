bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar maxDeltaEps =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaEps", 0.1);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

// ************************************************************************* //
