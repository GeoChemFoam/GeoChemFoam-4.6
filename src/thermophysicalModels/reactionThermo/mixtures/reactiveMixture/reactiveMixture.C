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

#include "reactiveMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactiveMixture::reactiveMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    basicMultiComponentMixture(thermoDict, specieNames, mesh, obj),
    kineticPhase_(thermoDict.subDict("kineticPhases").toc()),
    kineticPhaseReaction_(thermoDict.subDict("kineticPhaseReactions").toc()),
    Ke_(kineticPhase_.size()),
    k0_(kineticPhaseReaction_.size()),
    scoeff_(kineticPhaseReaction_.size()),
    kiSpeciesIndex_(kineticPhaseReaction_.size()),
    ki_(kineticPhaseReaction_.size()),
    Omega_(kineticPhase_.size()),
    R_(kineticPhaseReaction_.size()),
    component1Index_(NULL),
    component2Index_(NULL)
{
    component1Index_ = (int*)malloc((size_t)(kineticPhase_.size() * sizeof(int)));
    component2Index_ = (int*)malloc((size_t)(kineticPhase_.size() * sizeof(int)));
    
    dictionary kpDict = thermoDict.subDict("kineticPhases");
    forAll(kineticPhase_, i)
    {
        dictionary kpSubDict = kpDict.subDict(kineticPhase_[i]);
        speciesTable kpSpecies(kpSubDict.subDict("species").toc());
        forAll(species_,j)
        {
            if (species_[j]==kpSpecies[0])
            {
                component1Index_[i]=j;
                break;
            }
            if (j==species_.size())
            {
                Info<< kpSpecies[0] << "not find in solutionSpecies"
                << endl
                << abort(FatalError);
            }
        }
        
        forAll(species_,j)
        {
            if (species_[j]==kpSpecies[1])
            {
                component2Index_[i]=j;
                break;
            }
            if (j==species_.size())
            {
                Info<< kpSpecies[1] << "not find in solutionSpecies"
                << endl
                << abort(FatalError);
            }
        }

        
        Ke_.set
	(
		i,
		new dimensionedScalar(kpSubDict.lookup("Ke"))
	);
	
        Omega_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    kineticPhase_[i],
                    mesh.time().timeName(),
                    obj,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                Y_[component1Index_[i]]*Y_[component2Index_[i]]/Ke_[i]
            )
        );
    }
    
    dictionary kprDict = thermoDict.subDict("kineticPhaseReactions");
    forAll(kineticPhaseReaction_, j)
    {
        dictionary kprSubDict = kprDict.subDict(kineticPhaseReaction_[j]);
        k0_.set
	(
		j,
		new dimensionedScalar(kprSubDict.lookup("k0"))
	);
	
	scoeff_.set
	(
	       j,
	       new scalarList(species_.size(),0.0)
	);
	
	
	dictionary kprSpeciesDict = kprSubDict.subDict("species");
	speciesTable kprSpecies(kprSpeciesDict.toc());
	forAll(species_,i)
	{
	    if (kprSpeciesDict.isDict(species_[i]))
	    {
	        scoeff_[j][i] = readScalar(kprSpeciesDict.subDict(species_[i]).lookup("scoeff"));
	    }
	}

        
        kiSpeciesIndex_.set
        (
             j,
             new labelList(kprSpecies.size(),0)
        );
        
        forAll(kprSpecies,i)
        {
		forAll(species_,k)
		{
		    if (kprSpecies[i]==species_[k])
		    {
                        kiSpeciesIndex_[j][i]=k;
		    }
		}
	}
	
        ki_.set
	(
	      j,
	      new PtrList<dimensionedScalar>(kprSpecies.size())
	);
	
	forAll(kprSpecies,i)
        {
            ki_[j].set
            (
                i,
                new dimensionedScalar(kprSpeciesDict.subDict(kprSpecies[i]).lookup("ki"))
            );
        }	
	
	R_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "R_"+kineticPhaseReaction_[j],
                    mesh.time().timeName(),
                    obj,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                k0_[j]
            )
        );
    }

    correct();
}

Foam::reactiveMixture::reactiveMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& phaseName
)
:
    basicMultiComponentMixture(thermoDict, specieNames, mesh, obj, phaseName),
    kineticPhase_(thermoDict.subDict("kineticPhases").toc()),
    kineticPhaseReaction_(thermoDict.subDict("kineticPhaseReactions").toc()),
    Ke_(kineticPhase_.size()),
    k0_(kineticPhaseReaction_.size()),
    scoeff_(kineticPhaseReaction_.size()),
    kiSpeciesIndex_(kineticPhaseReaction_.size()),
    ki_(kineticPhaseReaction_.size()),
    Omega_(kineticPhase_.size()),
    R_(kineticPhaseReaction_.size()),
    component1Index_(NULL),
    component2Index_(NULL)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::reactiveMixture::~reactiveMixture
(
)
{
	free(component1Index_);
	free(component2Index_);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::reactiveMixture::correct()
{
    forAll(kineticPhase_, i)
    {
        Omega_[i] = Y_[component1Index_[i]]*Y_[component2Index_[i]]/Ke_[i];
    }
    
    forAll(kineticPhaseReaction_, j)
    {
        R_[j] =k0_[j];
        forAll(kiSpeciesIndex_[j],i)
        {
            R_[j] += ki_[j][i]*Y_[kiSpeciesIndex_[j][i]];
        }
    }
}

// ************************************************************************* //
