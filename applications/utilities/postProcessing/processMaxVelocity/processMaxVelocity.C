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
    processMaxVelocity 

Description
    calculate max velocity for each time-step

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "fvCFD.H"
#include "argList.H"
#include "primitivePatchInterpolation.H"
#include "timeSelector.H"


using namespace Foam;


int main(int argc, char *argv[])
{
	#   include "setRootCase.H"
	#   include "createTime.H"
	#   include "createNamedMesh.H"


	instantList timeList = timeSelector::select0(runTime, args);

	std::ofstream csvfile("maxVelocity.csv");
	csvfile << "time v\n";

	forAll(timeList, timeStep)
	{	
		runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;


			volVectorField U
			(
				IOobject
				(
					"U",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
					),
				mesh
				);


		scalar Um = max(mag(U)).value();  

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << Um << "\n";
		}
	}

	csvfile.close();

	return 0;
}


// ************************************************************************* //
