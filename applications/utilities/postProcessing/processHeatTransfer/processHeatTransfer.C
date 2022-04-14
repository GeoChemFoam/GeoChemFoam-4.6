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
    processheatTransfer 

Description
    calculate heat transfer coefficient between solid and fluid

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


	instantList timeList = timeSelector::select0(runTime, args);


        string filename = "HT.csv";
	std::ofstream csvfile(filename.c_str());
	csvfile << "time coeff\n";

        
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
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    )
                );

                scalar z1
                (
                    postProcessDict.lookupOrDefault<scalar>("z1",-1e30)
                );


                scalar z2
                (
                    postProcessDict.lookupOrDefault<scalar>("z2",1e30)
                );

                scalar y1
                (
                    postProcessDict.lookupOrDefault<scalar>("y1",-1e30)
                );


                scalar y2
                (
                    postProcessDict.lookupOrDefault<scalar>("y2",1e30)
                );


                scalar x1
                (
                    postProcessDict.lookupOrDefault<scalar>("x1",1e-30)
                );


                scalar x2
                (
                    postProcessDict.lookupOrDefault<scalar>("x2",1e30)
                );


                volScalarField clip
                (
                        IOobject
                        (
                            "clip",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh,
                        dimensionedScalar("clip",dimVolume,0.0)
                );

                volScalarField coordz=mesh.C().component(2);
                volScalarField coordy=mesh.C().component(1);
                volScalarField coordx=mesh.C().component(0);


               forAll(mesh.cells(),j)
               {
                    scalar zj=coordz[j];
                    scalar yj=coordy[j];
                    scalar xj=coordx[j];
                    if (zj>=z1 && zj<=z2 && yj>=y1 && yj<=y2 && xj>=x1 && xj<=x2)
                    {
                        clip[j]=mesh.V()[j];
                    }
                }

               Info<< "Reading transportProperties\n" << endl;

               IOdictionary transportProperties
               (
                   IOobject
                   (
                       "transportProperties",
                       runTime.constant(),
                       mesh,
                       IOobject::MUST_READ,
                       IOobject::NO_WRITE
                   )
               );


              Info<< "Reading diffusivity DT\n" << endl;

              dimensionedScalar DTf
              (
                  transportProperties.lookup("DTf")
              );

              dimensionedScalar DTs
              (
                  transportProperties.lookup("DTs")
              );

		volScalarField eps 
		(
			IOobject
			(
				"eps",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);

		volScalarField eps2 = 1-eps;

                volScalarField T
                (
                        IOobject
                        (
                                "T",
                                runTime.timeName(),
                                mesh,
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE
                        ),
                        mesh
                );


                scalar Tf = T.weightedAverage(eps*clip).value();

                Info << "Tf:" << Tf << endl;

                scalar Ts = T.weightedAverage(eps2*clip).value();

                Info << "Ts:" << Ts << endl;


                surfaceScalarField D = DTf*DTs/fvc::interpolate(eps*DTs + (1-eps)*DTf);

            
                surfaceScalarField phiD = D*fvc::snGrad(T)*mesh.magSf();

                surfaceScalarField epsf = fvc::interpolate(eps);
 
                volScalarField F = fvc::div(phiD*epsf)-eps*fvc::div(phiD);

                scalar HT = F.weightedAverage(clip).value();

                Info << "HT:" << HT << endl;

                volScalarField gEps = mag(fvc::grad(eps));

                scalar A = gEps.weightedAverage(clip).value();

                Info << "A:" << A << endl;
 

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << -HT/A/(Ts-Tf)  << "\n";
		}

        }

        csvfile.close();

	return 0;
}


// ************************************************************************* //
