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

	instantList timeList = timeSelector::select0(runTime, args);


	forAll(timeList, timeStep)
	{	
		runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;

		#   include "createNamedMesh.H"

                volVectorField cellCenters
                (
                        IOobject
                        (
                                "cellCenters",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh.C()
                );

                cellCenters.write();
	}

	return 0;
}


// ************************************************************************* //
