    // Solve the momentum equation

    tmp<fvVectorMatrix> HUEqn
    (
      (1/eps/eps)*fvm::div(phi, U)
      - 1/eps*fvm::laplacian(nu, U)
      + fvm::Sp(nu*Kinv,U)
    );

    // Get under-relaxation factor
    const scalar UUrf = mesh.solutionDict().equationRelaxationFactor(U.name());

    // Momentum solution
    if (steadyState.momentumPredictor())
    {
        solve
        (
            relax(HUEqn(), UUrf)
         ==
            -fvc::grad(p)
        );
    }
