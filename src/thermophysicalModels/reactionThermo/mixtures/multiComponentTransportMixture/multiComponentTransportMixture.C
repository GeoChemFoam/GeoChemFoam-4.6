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

#include "multiComponentTransportMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::multiComponentTransportMixture<MixtureType>::multiComponentTransportMixture
(
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    IOdictionary
    (
        IOobject
        (
            "thermoPhysicalProperties",
            mesh.time().constant(),
            obj,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    MixtureType(*this, this->subDict("solutionSpecies").toc(), mesh, obj),
	DY_(this->subDict("solutionSpecies").toc().size())
{
    wordList specieNames(this->subDict("solutionSpecies").toc());
	Info << "Read species diffusion coefficients\n" << endl;
	forAll(specieNames, i)
	{
	    const dictionary& solutionSpeciesDict = this->subDict("solutionSpecies");
		const dictionary& subdict = solutionSpeciesDict.subDict(specieNames[i]);
		DY_.set
		(
			i,
			new dimensionedScalar(subdict.lookup("D"))
		);
	}
}

template<class MixtureType>
Foam::multiComponentTransportMixture<MixtureType>::multiComponentTransportMixture
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "thermoPhysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    MixtureType(*this, this->subDict("solutionSpecies").toc(), mesh, mesh),
	DY_(this->subDict("solutionSpecies").toc().size())
{
    wordList specieNames(this->subDict("solutionSpecies").toc());
	Info << "Read species diffusion coefficients\n" << endl;
	forAll(specieNames, i)
	{
	    const dictionary& solutionSpeciesDict = this->subDict("solutionSpecies");
		const dictionary& subdict = solutionSpeciesDict.subDict(specieNames[i]);
		DY_.set
		(
			i,
			new dimensionedScalar(subdict.lookup("D"))
		);
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
namespace Foam
{
    template class Foam::multiComponentTransportMixture<Foam::basicMultiComponentMixture>;
    template class Foam::multiComponentTransportMixture<Foam::solutionSurfaceMultiComponentMixture>;
    template class Foam::multiComponentTransportMixture<Foam::phreeqcMixture>;
    template class Foam::multiComponentTransportMixture<Foam::reactiveMixture>;
}

