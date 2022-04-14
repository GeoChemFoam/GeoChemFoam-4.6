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
    reactiveTransportDBSFoam

Description
    Solves reactive transport equation with microcontinuum DBS for multi-species flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "multiComponentTransportMixture.H"
#include "steadyStateControl.H"
#include "reactiveMixture.H"

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
#   include "deltaEpsMax.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "deltaEpsMax.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

#       include "setDeltaT.H"
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

#       include "volContinuity.H"

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);
        
#       include "epsEqn.H"

        //permeability
        Kinv = kf*pow(1-eps,2)/pow(eps,3);

        volVectorField gradEps = fvc::grad(eps);
        surfaceVectorField gradEpsf = fvc::interpolate(gradEps);
        surfaceVectorField nEpsv = -gradEpsf/(mag(gradEpsf) + deltaN);
        nEpsf = nEpsv & mesh.Sf();

        volScalarField a = mag(fvc::grad(eps));

        scalar lambda = psiCoeff;

        if (VoS=="VoS-psi")
        {
            scalar As = a.weightedAverage(mesh.V()).value();

            a = a*(1-eps)*(1e-3+eps);

            if (adaptPsiCoeff) lambda = As/a.weightedAverage(mesh.V()).value();

            a = lambda*a;


            Info << "psiCoeff=" << lambda << endl;
       }


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
