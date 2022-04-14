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
    processPoroPerm

Description
    calculate porosity and permeability from flow field for each time-step

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

	std::ofstream csvfile("poroPerm.csv");
	csvfile << "time poro perm\n";

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

		IOdictionary transportProperties 
		(
			IOobject
			(
				"transportProperties",
				"constant",
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);

		dimensionedScalar nu
		(
		    transportProperties.lookup("nu")
		);


                word velocity
                (
                    postProcessDict.lookup("velocity")
                );
 
                word pressureDrop 
                (
                    postProcessDict.lookup("pressureDrop")
                );


		label direction
		(
		    readLabel(postProcessDict.lookup("direction"))
		);

                scalar z1
                (
                    readScalar(postProcessDict.lookup("z1"))
                );


                scalar z2
                (
                    readScalar(postProcessDict.lookup("z2"))
                );

                scalar y1
                (
                    readScalar(postProcessDict.lookup("y1"))
                );


                scalar y2
                (
                    readScalar(postProcessDict.lookup("y2"))
                );


                scalar x1
                (
                    readScalar(postProcessDict.lookup("x1"))
                );


                scalar x2
                (
                   readScalar(postProcessDict.lookup("x2"))
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

                scalar vol = gSum(clip);


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

                volScalarField p
                (
                        IOobject
                        (
                                "p",
                                runTime.timeName(),
                                mesh,
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh
                );

                scalar L=0.0;

                if (direction==0) L=x2-x1;
                if (direction==1) L=y2-y1;
                if (direction==2) L=z2-z1;

                volScalarField Ud = U.component(direction);
                scalar Uavg =  Ud.weightedAverage(clip).value();
                scalar UD=Uavg*vol/(x2-x1)/(y2-y1)/(z2-z1);


               if (velocity=="inlet")
               {
                   word inletName
                   (
                       postProcessDict.lookup("inletName")
                   );

                   label inletPatchID  = mesh.boundaryMesh().findPatchID(inletName);

                   fvPatchVectorField& Uf = U.boundaryField()[inletPatchID];
                   scalarField Ufd = Uf.component(direction);

                   fvPatchScalarField& pInlet  = p.boundaryField()[inletPatchID];
                   scalarField magSfInlet  = pInlet.patch().magSf();

                   scalar Q  = gSum(Ufd*magSfInlet);
                   scalar A  = (x2-x1)*(y2-y1)*(z2-z1)/L;
                   UD = Q/A;

               }
               else if (velocity=="volumeAveraged")
               {
                  //UD=Uavg*vol/(x2-x1)/(y2-y1)/(z2-z1);
               }
               else
               { 
                   Info<< "velocity should be equal to inlet or volumeAveraged"
                       << endl
                       << abort(FatalError);
               }

                volScalarField gradRate= fvc::grad(p) & U;
                scalar gradP = gradRate.weightedAverage(clip).value();

                scalar dpdz= -gradP/(Uavg+1e-30);


                if (pressureDrop=="direct")
                {
                    word inletName
                    (
                        postProcessDict.lookup("inletName")
                    );

                    word outletName
                    (
                        postProcessDict.lookup("outletName")
                    );

                    label inletPatchID  = mesh.boundaryMesh().findPatchID(inletName);
                    label outletPatchID = mesh.boundaryMesh().findPatchID(outletName);

                    fvPatchScalarField& pInlet  = p.boundaryField()[inletPatchID];
                    fvPatchScalarField& pOutlet = p.boundaryField()[outletPatchID];
                    scalarField magSfInlet  = pInlet.patch().magSf();
                    scalarField magSfOutlet = pOutlet.patch().magSf();

                    dpdz = mag(gSum(pInlet*magSfInlet)/gSum(magSfInlet)-gSum(pOutlet*magSfOutlet)/gSum(magSfOutlet))/L+1e-13;

                }
                else if (pressureDrop=="velocityWeightedAverageGradient")
                {
                    //dpdz= -gradP/(Uavg+1e-30);
                }
                else if (pressureDrop=="velocityWeightedAverageViscousForce")
                {
                    volScalarField viscPRate=(nu/eps)*(fvc::laplacian(U) & U);

                    scalar viscP = viscPRate.weightedAverage(clip).value();

                    dpdz=-viscP/(Uavg+1e-30);
                }
               else
               {
                   Info<< "pressure drop should be equal to direct or velocityWeightedAverageGradient or velocityWeightedAverageViscousForce"
                       << endl
                       << abort(FatalError);
               }



		scalar poro = eps.weightedAverage(clip).value()*vol/(x2-x1)/(y2-y1)/(z2-z1);

                scalar perm = UD*nu.value()/(dpdz+1e-30);
                //scalar perm = -Uavg*vol/(x2-x1)/(y2-y1)/(z2-z1)*nu.value()/(dpdz+1e-30); 


		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << poro << " " << perm << "\n";
		}
	}

	csvfile.close();

	return 0;
}


// ************************************************************************* //
