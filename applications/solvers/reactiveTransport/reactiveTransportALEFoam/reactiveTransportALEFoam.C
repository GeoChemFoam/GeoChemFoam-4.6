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
    reactiveTransportALEFoam

Description
    Solves reactive transport equation with ALE mesh for multi-species flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "multiComponentTransportMixture.H"
#include "reactiveMixture.H"
#include "steadyStateControl.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    simpleControl simple(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

#   include "createTimeControls.H"
#   include "meshCourantNo.H"
#   include "setInitialDeltaT.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "meshCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.update();

    	steadyStateControl steadyState(mesh);
	while (steadyState.loop())
	{
                // Pressure-velocity SIMPLE corrector
		{
#               include "UEqn.H"
#               include "pEqn.H"
		}
		

	        turbulence->correct();

               // Concentration solver
               {
#               include "YiEqn.H"
               }

        }

                if (!steadyState.criteriaSatisfied())
                {
                    U == U.oldTime();
                    p == p.oldTime();
                    forAll(solutionSpecies,i)
                    {
                        volScalarField Yi = speciesMixture.Y(i);
                        Yi == Yi.oldTime(); 
                    }
                     //runTime.writeNow();
                     //return 0;
                }

                if (runTime.timeIndex() > remeshFrequency-1 )
                {
                     runTime.writeNow();
                     return 0;
                }


		if (runTime.timeIndex() % checkFrequency == 0)
		{
			if (mesh.checkMesh(false))
			{
				runTime.writeNow();
				return 0;
			}
		}

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    	runTime.write();
    }

    runTime.writeNow();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
