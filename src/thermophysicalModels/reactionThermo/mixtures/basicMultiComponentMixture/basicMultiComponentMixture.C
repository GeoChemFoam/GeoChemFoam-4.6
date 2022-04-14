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

#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMultiComponentMixture::basicMultiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    species_(specieNames),
    Y_(species_.size())
{
    forAll(species_, i)
    {
        IOobject header
        (
            species_[i],
            mesh.time().timeName(),
            obj,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        obj,
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
                    obj,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        obj,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault
                )
            );
        }
    }

    // Do not enforce constraint of sum of mass fractions to equal 1 here
    // - not applicable to all models
}

Foam::basicMultiComponentMixture::basicMultiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& phaseName
)
:
    species_(specieNames),
    Y_(species_.size())
{
    forAll(species_, i)
    {
        Y_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    species_[i]+"_"+phaseName,
                    mesh.time().timeName(),
                    obj,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
		        dimensionedScalar(species_[i]+"_"+phaseName, dimMoles/dimVolume, 0.0)
            )
        );
    }
}


// ************************************************************************* //
