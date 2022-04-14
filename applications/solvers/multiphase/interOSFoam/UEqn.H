    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho_vd*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_vd, Uvd)
      + fvm::div(rhoPhi, Uvd, "div(rho*phi,U)")
      - fvm::laplacian(muEff, Uvd)
      - (fvc::grad(Uvd) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  //  fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho_vd)
                  - fvc::snGrad(pvd)
                )*mesh.magSf()
            )
        );
    }
