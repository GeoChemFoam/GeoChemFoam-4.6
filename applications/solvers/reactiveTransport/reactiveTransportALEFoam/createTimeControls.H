bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar maxMeshCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxMeshCo", 1);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

// ************************************************************************* //
