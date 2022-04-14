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
    processPhaseConcentration 

Description
    calculate average concentration of each species in each phase for each time-step

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "fvCFD.H"
#include "argList.H"
#include "primitivePatchInterpolation.H"
#include "timeSelector.H"
#include "twoPhaseMixture.H"
#include "speciesTable.H"


using namespace Foam;


int main(int argc, char *argv[])
{
	#   include "setRootCase.H"
	#   include "createTime.H"
	#   include "createNamedMesh.H"


	instantList timeList = timeSelector::select0(runTime, args);

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
        IOdictionary thermoPhysicalProperties
        (
                IOobject
                (
                       "thermoPhysicalProperties",
                       "constant",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                )
        );

	const dictionary& solutionSpeciesDict = thermoPhysicalProperties.subDict("solutionSpecies");
	speciesTable solutionSpecies
	(
		solutionSpeciesDict.toc()
	);


	PtrList<dimensionedScalar> HY(solutionSpecies.size());

	forAll(solutionSpecies, i)
	{
		const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		HY.set
		(
			i,
			new dimensionedScalar(subdict.lookup("H"))
		);
	}


        forAll(solutionSpecies,i)
        {
            string filename = solutionSpecies[i]+"_phaseConc.csv";
	    std::ofstream csvfile(filename.c_str());
	    csvfile << "time C1 C2\n";

        
	    forAll(timeList, timeStep)
	    {	
	    	runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;

		volScalarField alpha1 
		(
			IOobject
			(
				"alpha1",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);

		volScalarField alpha2 = 1-alpha1;

              	PtrList<volScalarField> Y(solutionSpecies.size());

                IOobject header
		(
			solutionSpecies[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (header.headerOk())
		{
			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			volScalarField Ydefault
			(
				IOobject
				(
					"Ydefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					Ydefault
				)
			);
		}

		volScalarField& Yi = Y[i];
		dimensionedScalar HYi = HY[i];
		volScalarField Y1i = alpha1*Yi/(alpha1+HYi*(1-alpha1));
		volScalarField Y2i = alpha2*HYi*Yi/(alpha1+HYi*(1-alpha1));

                scalarField& vol = clip.internalField();
            
		scalar T1 = gSum(Y1i*vol)/(gSum(alpha1*vol)+1e-30);
		scalar T2 = gSum(Y2i*vol)/(gSum(alpha2*vol)+1e-30);

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << T1 << " " << T2 << "\n";
		}
	    }

	    csvfile.close();
        }

	return 0;
}


// ************************************************************************* //
