{
    volScalarField rUA("rUA",1.0/UEqn.A());
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    Uvd = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(Uvd) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho_vd, Uvd, phivd)
    );

    adjustPhi(phiU, Uvd, pvd);

    phivd = phiU +
        (
          //  fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho_vd)
        )*rUAf*mesh.magSf();


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pvd) == fvc::div(phivd)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phivd -= pdEqn.flux();
        }
    }

    Uvd += rUA*fvc::reconstruct((phivd - phiU)/rUAf);
    Uvd.correctBoundaryConditions();
}
