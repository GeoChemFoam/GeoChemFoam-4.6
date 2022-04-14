{
    volScalarField rUA("rUA",1.0/URelaxEqn.A());
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    Ucr = rUA*URelaxEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(Ucr) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho_cr, Ucr, phicr)
    );

    adjustPhi(phiU, Ucr, pcr);

    phicr = phiU +
        rUAf*interface.phic();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pcr) == fvc::div(phicr)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );


        if (pimple.finalNonOrthogonalIter())
        {
            phicr -= pdEqn.flux();
        }
    }

    Ucr += rUA*fvc::reconstruct((phicr - phiU)/rUAf);
    Ucr.correctBoundaryConditions();
}
