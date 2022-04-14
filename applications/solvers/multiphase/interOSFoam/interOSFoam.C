/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the Carbonate Reservoir Group at Heriot-Watt University. Please visit our
    website for more information <https://carbonates.hw.ac.uk>.

Application
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e.  laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//define monitor time indexes
#define ifMonitor  if( runTime.timeIndex()%10== 0 )

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pimpleControl pimple(mesh);

#   include "readGravitationalAcceleration.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "createTimeControls.H"
#   include "correctPhivd.H"
#   include "correctPhicr.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity corrector
        while (pimple.loop())
        {
            twoPhaseProperties.correct();

#           include "alphaEqnSubCycle.H"

            U = 0*U;
            phi = 0*phi;

            // Initialize alpha for capillary relaxation step
            alpha1.oldTime() = alpha1;
                        int nRelax = ::floor(runTime.deltaTValue()/deltaTCR-1e-9)+1;
            runTime.setDeltaT(runTime.deltaTValue()/nRelax,false);

            nSteps++;
            nRelaxTot += nRelax;


            for (int i=0;i<nRelax;i++)
            { 
                nRelaxSteps++;


                Info << "Time=" << runTime.timeName() << ", relaxation iteration:" << i+1 << " out of " << nRelax << endl;

                twoPhaseProperties.correct();

#               include "alphaRelaxEqnSubCycle.H"

                rhoPhi = rhoPhivd + rhoPhicr;

#               include "URelaxEqn.H"

                // --- PISO loop
                while (pimple.correct())
                {
#                   include "pRelaxEqn.H"
                }

#               include "continuityErrsRelax.H"

                Info << "\n  Ucrmax = " << max(mag(Ucr)).value() << " m/s  " ;


                alpha1.oldTime() = alpha1;
                Ucr.oldTime() = Ucr;
                rho_cr.oldTime() = rho_cr;
                phicr.oldTime() = phicr;

                U += Ucr/nRelax;
                phi += phicr/nRelax;


                if (i>=relaxMin-1 && i> 0.9*(nIterLast-1) && interface.eqnResidual() < residual) 
                {
                    Info << "\n         Time=" << runTime.timeName() << ", capillary relaxation finished in " << 1+i << " steps out of " << nRelax << endl;
                    nIterLast = i;
                    break;
                }
                if (i==nRelax-1)
                {
                    Info << "\n         Time=" << runTime.timeName() << ", capillary relaxation finished in " << 1+i << " steps out of " << nRelax << endl;
                    nIterLast = i;
                }
            }

            Info << "\n         Time=" << runTime.timeName() << ", end of capillary relaxation" << endl; 
            Info << "\n         Average number of relaxation steps="<< nRelaxSteps/nSteps << endl;
            Info << "\n         Percent of relaxation steps performed="<< nRelaxSteps/nRelaxTot*100 << "%" << endl;
            Info << "\n" << endl;

            runTime.setDeltaT(nRelax*runTime.deltaTValue(),false);

            rho_vd = rho_cr;

#           include "UEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
#               include "pEqn.H"
            }

#           include "continuityErrs.H"


            pd = pvd + pcr;

            p = pd + pc + rho_vd*gh;

            if (pd.needReference())
            {
                p += dimensionedScalar
                (
                    "p",
                    p.dimensions(),
                    pRefValue - getRefCellValue(pd, pdRefCell)
                );
            }

            U += Uvd;
            phi += phivd;

            turbulence->correct();
        }

        runTime.write();

        //monitor average and max velocity
        ifMonitor
        {
            Info << "\n         Umax = " << max(mag(U)).value() << " m/s  "
            << "Uavg = " << mag(average(U)).value() << " m/s";
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
