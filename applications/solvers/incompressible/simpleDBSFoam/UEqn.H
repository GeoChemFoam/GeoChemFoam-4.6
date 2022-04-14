    // Solve the momentum equation


    tmp<fvVectorMatrix> HUEqn
    (
        (1/eps/eps)*fvm::div(phi, U)
      + (1/eps)*turbulence->divDevReff()
      + fvm::Sp(nu*Kinv,U)
      - momentumSource 
    );

    // Get under-relaxation factor
    const scalar UUrf = mesh.solutionDict().equationRelaxationFactor(U.name());

    // Momentum solution
    solve
    (
        relax(HUEqn(), UUrf)
     ==
        -fvc::grad(p)
    );
