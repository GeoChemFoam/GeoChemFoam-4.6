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

\*---------------------------------------------------------------------------*/

#include "phreeqcMixture.H"
#include "RM_interface_C.h"
#include "fvcAverage.H"
#include "surfaceInterpolate.H"
#include "reactingWallFvPatchScalarField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::phreeqcMixture::initialise()
{
	//get number of cells
	int ncells = mesh_.cells().size();
	id_ = RM_Create(ncells, 1);

	//set porosity=1.0
	double* porosity = new double[ncells];
	for (int i = 0; i < ncells; i++) porosity[i] = 1.0;
	RM_SetPorosity(id_, porosity);
	delete[] porosity;

	
	//concentration=mol/L
	RM_SetUnitsSolution(id_, 2);

	//Save species for multi-species transport
	RM_SetSpeciesSaveOn(id_, true);

	//load phreeqc database
	RM_LoadDatabase(id_, "constant/GeoChem.dat");

	//load reaction
	RM_RunFile(id_, 1, 1, 0, "constant/phreeqcReactions");

	//solution component list
	std::ostringstream oss;


	//init surface and solution
	forAll(mesh_.cells(),celli)
	{
		oss << "COPY solution 0 " << celli  << "\n";
		oss << "END" << "\n";
    }

	//equilibrate surface with solution
	forAll(mesh_.cells(),celli)
	{
		oss << "SURFACE " << celli << "\n";
		oss << "-equilibrate with solution " << celli << "\n";
		forAll(surfaceMasters_, j)
		{
			double area = 0.0;
			double mole = 0.0;
	        const volScalarField::GeometricBoundaryField& Surfbf = Surf_.boundaryField();
			forAll(Surfbf,patchi)
			{
				if (Surfbf[patchi].type() == "reactingWall")
				{
					const reactingWallFvPatchScalarField& Surfcap = refCast<const reactingWallFvPatchScalarField>(Surfbf[patchi]);
					const labelList& cellOwner = Surfcap.patch().faceCells();
					const surfaceScalarField& magSf = mesh_.magSf();
					const wordList& masters = Surfcap.get_surface_masters();
					const scalarList& density   = Surfcap.get_density();//mol/m^2
					forAll(masters,i)
					{
						if (masters[i]==surfaceMasters_[j])
						{
							forAll(Surfbf[patchi],facei)
							{
								if (cellOwner[facei]==celli)
								{
									mole+=density[i]*magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]];//mol/L
									area+=magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;//m^2/L
								}
							}
						}
					}
				}
			}	
            if (area>0) oss << surfaceMasters_[j] << "  " << mole << " " << area  << " 1" << "\n";
            else oss << surfaceMasters_[j] << "  " << "0 1 1" << "\n";
		}
		oss << "END" << "\n";
	}

	//run phreeqc keywords
	RM_RunString(id_, 1, 1, 0, oss.str().c_str());

	//init Phreeqc worker module
	int* ic1 = (int *)malloc((size_t)(7 * ncells * sizeof(int)));
	forAll(mesh_.cells(),celli)
	{
		ic1[celli] = celli;               // Solution  i
		ic1[ncells + celli] = -1;      // Equilibrium phases none
		ic1[2 * ncells + celli] = -1;       // Exchange none
		ic1[3 * ncells + celli] = celli;      // Surface i
		ic1[4 * ncells + celli] = -1;      // Gas phase none
		ic1[5 * ncells + celli] = -1;      // Solid solutions none
		ic1[6 * ncells + celli] = -1;      // Kinetics none
	}
	RM_InitialPhreeqc2Module(id_, ic1, 0, 0);
	free(ic1);

	//Run Phreeqc to init concentration
	RM_RunCells(id_);

	//find components
	RM_FindComponents(id_);
	
	//get number of solution species
	int nsol = RM_GetSpeciesCount(id_);

	//display number of silution species
	Info << "nsol:" << nsol << "\n";

	//set solution speciesvector
	char** components = (char **)malloc((size_t)(nsol * sizeof(char *)));

	//get solution species name
	for (int i = 0; i < nsol; i++)
	{
		components[i] = (char *)malloc((size_t)(20 * sizeof(char *)));
		RM_GetSpeciesName(id_, i, components[i], 20);
	}

	//get number of surface species
	int nsurf = RM_GetSurfaceSpeciesCount(id_);

	//display number of surface species
	Info << "nsurf:" << nsurf << "\n";

	//set solution speciesvector
	char** surfComponents = (char **)malloc((size_t)(nsurf * sizeof(char *)));

	//get solution species name
	for (int i = 0; i < nsurf; i++)
	{
		surfComponents[i] = (char *)malloc((size_t)(20 * sizeof(char *)));
		RM_GetSurfaceSpeciesName(id_, i, surfComponents[i], 20);
	}
	
	//save component index map
	componentSolutionIndex_ = (int*)malloc((size_t)(species_.size() * sizeof(int)));
	componentSurfaceIndex_  = (int*)malloc((size_t)(surfaceSpecies_.size() * sizeof(int)));
	forAll(species_, i)
	{
		for (int j = 0; j < nsol; j++)
		{
			std::string component = components[j];
			if (component == species_[i])
			{
				componentSolutionIndex_[i] = j;
			}
		}
	}

	forAll(surfaceSpecies_, i)
	{
		for (int j = 0; j < nsurf; j++)
		{
			std::string component = surfComponents[j];
			if (component == surfaceSpecies_[i])
			{
				componentSurfaceIndex_[i] = j;
			}
		}
	}

	//concentration, adsorption and surface potential
	if (nsol>0) concentration_      = (double *)malloc((size_t)(nsol * ncells * sizeof(double)));
	if (nsurf>0) 
	{
		surfConcentration_  = (double *)malloc((size_t)(nsurf * ncells * sizeof(double)));
		surfArea_           = (double *)malloc((size_t)(ncells * sizeof(double)));
		surfPotential_      = (double *)malloc((size_t)(ncells * sizeof(double)));
	}
	
	RM_GetSpeciesConcentrations(id_, concentration_);

	RM_GetSurfaceSpeciesConcentrations(id_, surfConcentration_);
	RM_GetSurfaceArea(id_,"Surf",surfArea_);
	RM_GetSurfacePotential(id_,"Surf",surfPotential_);

	//set water saturation for Phreeqc module
	saturation_ = (double *)malloc((size_t)(ncells * sizeof(double)));
    forAll(mesh_.cells(),celli)
    {
        saturation_[celli]=1.0;
    }
    //memset(saturation_,1.0,ncells * sizeof(double));
	RM_SetSaturation(id_, saturation_);

	//get concentration from OpenFOAM
	forAll(species_, i)
	{
		volScalarField& Yi = Y_[i];
		forAll(mesh_.cells(), celli)
		{
			concentration_[componentSolutionIndex_[i] * ncells + celli] = Yi[celli];
		}
	}

	//get surface concentration from OpenFoam
	forAll(surfaceSpecies_, i)
	{
		forAll(mesh_.cells(), celli)
		{
			surfConcentration_[componentSurfaceIndex_[i] * ncells + celli] = 0.0;
		}
		volScalarField& Yi = surfaceMixture_.Y(i);
	    forAll(Yi.boundaryField(), patchi)
		{
			if (Surf_.boundaryField()[patchi].type()=="reactingWall")
			{
				const labelList& cellOwner = Yi.boundaryField()[patchi].patch().faceCells();
				const scalarField& Yfaces = Yi.boundaryField()[patchi];
				const surfaceScalarField& magSf = mesh_.magSf();
				forAll(Yi.boundaryField()[patchi], facei)
				{
					if (surfArea_[cellOwner[facei]]>0)
					{
						surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]] += Yfaces[facei] / surfArea_[cellOwner[facei]] * magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]];
					}
				}
			}
		}
	}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phreeqcMixture::phreeqcMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    solutionSurfaceMultiComponentMixture(thermoDict, specieNames, mesh, obj),
    mesh_(mesh),
    Surf_
    (
        IOobject
        (
            "Surf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("Surf", dimMoles/dimArea, 0.0)
    ),
    I_
    (
        IOobject
        (
            "I",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	    dimensionedScalar("I", dimMoles/dimVolume, 0.0)
    ),
    psi_
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
	    dimensionedScalar("psi", dimensionSet(1,2,-3,0,0,-1,0), 0.0)
    ),
    saturation_(NULL),
    componentSolutionIndex_(NULL),
    componentSurfaceIndex_(NULL),
    concentration_(NULL),
    surfConcentration_(NULL),
    surfArea_(NULL),
    surfPotential_(NULL)
{
    initialise();
}

Foam::phreeqcMixture::phreeqcMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& phaseName
)
:
    solutionSurfaceMultiComponentMixture(thermoDict, specieNames, mesh, obj, phaseName),
    mesh_(mesh),
    Surf_
    (
        IOobject
        (
            "Surf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("Surf", dimMoles/dimArea, 0.0)
    ),
    I_
    (
        IOobject
        (
            "I",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	    dimensionedScalar("I", dimMoles/dimVolume, 0.0)
    ),
    psi_
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
	    dimensionedScalar("psi", dimensionSet(1,2,-3,0,0,-1,0), 0.0)
    ),
    saturation_(NULL),
    componentSolutionIndex_(NULL),
    componentSurfaceIndex_(NULL),
    concentration_(NULL),
    surfConcentration_(NULL),
    surfArea_(NULL),
    surfPotential_(NULL)
{
    initialise();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::phreeqcMixture::~phreeqcMixture
(
)
{
	free(saturation_);
	free(componentSolutionIndex_);
	free(componentSurfaceIndex_);
	free(concentration_);
	free(surfConcentration_);
	free(surfArea_);
	free(surfPotential_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::phreeqcMixture::correct()
{
	//get number of cells
	int ncells = mesh_.cells().size();

	//get concentration after transport
	forAll(species_, i)
	{
		volScalarField& Yi = Y_[i];
		forAll(mesh_.cells(), celli)
		{
			concentration_[componentSolutionIndex_[i] * ncells + celli] = Yi[celli];//mol/L
		}
	}


	//get surface concentration after transport
	forAll(surfaceSpecies_, i)
	{
		forAll(mesh_.cells(), celli)
		{
			surfConcentration_[componentSurfaceIndex_[i] * ncells + celli] = 0.0;
		}
		volScalarField& Yi = surfaceMixture_.Y(i);
	    forAll(Yi.boundaryField(), patchi)
		{
			if (Surf_.boundaryField()[patchi].type()=="reactingWall")
			{
				const labelList& cellOwner = Yi.boundaryField()[patchi].patch().faceCells();
				const scalarField& Yfaces = Yi.boundaryField()[patchi];
				const surfaceScalarField& magSf = mesh_.magSf();
				forAll(Yi.boundaryField()[patchi], facei)
				{
					if (surfArea_[cellOwner[facei]]>0)
					{
						surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]] += Yfaces[facei] / surfArea_[cellOwner[facei]] * magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]];
					}
				}
			}
		}
	}

	//set concentration for Phreeqc module
	RM_SpeciesConcentrations2Module(id_, concentration_);
	RM_SurfaceSpeciesConcentrations2Module(id_, surfConcentration_);

	//Run phreeqc
	RM_RunCells(id_);

	//get concentration after reactions
	RM_GetSpeciesConcentrations(id_, concentration_);
	RM_GetSurfaceSpeciesConcentrations(id_, surfConcentration_);

	//get surface potential
	RM_GetSurfacePotential(id_,"Surf",surfPotential_);

	//get ionic strength
	RM_GetSolutionIonicStrength(id_, &I_[0]);

	//get reaction rate and new vector composition
	forAll(species_, i)
	{
		volScalarField& Yi = Y_[i];
		forAll(mesh_.cells(), celli)
		{
			if (saturation_[celli]>1e-3)
			{
				Yi[celli]  = concentration_[componentSolutionIndex_[i] * ncells + celli];//mol/m3
			}
		}
	}

	//get surface concentration
	forAll(surfaceSpecies_, i)
	{
		volScalarField& Yi = surfaceMixture_.Y(i);
		forAll(Yi.boundaryField(), patchi)
		{
			if (Surf_.boundaryField()[patchi].type()=="reactingWall")
			{
				const labelList& cellOwner = Yi.boundaryField()[patchi].patch().faceCells();
				scalarField& Yfaces   = Yi.boundaryField()[patchi];
				scalarField& psifaces = psi_.boundaryField()[patchi];
				forAll(Yi.boundaryField()[patchi], facei)
				{
					if (saturation_[cellOwner[facei]]>1e-3)
					{
						Yfaces[facei] = surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]]/1000;//kmol/m2
						psifaces[facei]   = surfPotential_[cellOwner[facei]];
					}
				}
			}
		}
	}
}

void Foam::phreeqcMixture::setSaturation(const volScalarField& alpha)
{
	forAll(mesh_.cells(), celli)
	{
		if (alpha[celli] >1e-3) saturation_[celli] = alpha[celli];
		else  saturation_[celli]=0;
	}
	RM_SetSaturation(id_, saturation_);
}
// ************************************************************************* //
