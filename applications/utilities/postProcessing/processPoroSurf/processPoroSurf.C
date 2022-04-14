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

	std::ofstream csvfile("poroSurf.csv");
	csvfile << "time poro a\n";

	forAll(timeList, timeStep)
	{
		runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;
	
		#   include "createNamedMesh.H"

		IOdictionary postProcessDict 
		(
			IOobject
			(
				"postProcessDict",
				"system",
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);



		word solidName
		(
			postProcessDict.lookup("solidName")
		);


		scalar dxmin
		(
			readScalar(postProcessDict.lookup("x1"))
		);

		scalar dxmax
		(
			readScalar(postProcessDict.lookup("x2"))
		);

		scalar dymin
		(
			readScalar(postProcessDict.lookup("y1"))
		);

		scalar dymax
		(
			readScalar(postProcessDict.lookup("y2"))
		);

		scalar dzmin
		(
			readScalar(postProcessDict.lookup("z1"))
		);

		scalar dzmax
		(
			readScalar(postProcessDict.lookup("z2"))
		);

		label solidPatchID  = mesh.boundaryMesh().findPatchID(solidName);

		volScalarField eps 
		(
			IOobject
			(
				"eps",
				runTime.timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			mesh,
		    dimensionedScalar("eps",dimless,1.0)
		);

               scalar a = 0.0;

                if (solidPatchID>-1)
                {
		    scalarField magSfSolid  = mesh.boundary()[solidPatchID].magSf();

		    a = gSum(magSfSolid);
                }

		scalar poro = eps.weightedAverage(mesh.V()).value()*gSum(mesh.V())/(dxmax-dxmin)/(dymax-dymin)/(dzmax-dzmin);

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << poro << " " << a << "\n";
		}
	}

	csvfile.close();

	return 0;
}


// ************************************************************************* //
