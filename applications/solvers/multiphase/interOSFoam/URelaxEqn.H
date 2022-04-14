    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho_cr*turbulence->nut())
    );


    fvVectorMatrix URelaxEqn
    (
        fvm::ddt(rho_cr, Ucr)
      + fvm::div(rhoPhi, Ucr,"div(rho*phi,U)")
      - fvm::laplacian(muEff, Ucr)
      - (fvc::grad(Ucr) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    URelaxEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            URelaxEqn
         ==
            fvc::reconstruct
            (
                (
                  //  fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - fvc::snGrad(pcr)
                ) * mesh.magSf()
				//capillary flux
				+ interface.phic()
            )
        );
    }
