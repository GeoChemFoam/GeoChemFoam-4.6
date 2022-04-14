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

Class
    Foam::inertMultiComponentMixture

Description
    Foam::inertMultiComponentMixture

SourceFiles
    inertMultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef inertMultiComponentMixture_H
#define inertMultiComponentMixture_H

#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class inertMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

class inertMultiComponentMixture
:
    public basicMultiComponentMixture
{

protected:
    // Private member functions

        //- Construct as copy (not implemented)
        inertMultiComponentMixture(const inertMultiComponentMixture&);


public:


    // Constructors

        //- Construct from mesh and object registry
        inertMultiComponentMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&
        );

        //- Construct from mesh and object registry and phase name
        inertMultiComponentMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&,
            const word&
        );

        //- Construct from mesh
        inertMultiComponentMixture
        (
            const fvMesh&
        );


    // Member functions
    void correct()
	{
        //nothing to do
	}

    void setSaturation(const volScalarField& alpha)
	{
        //nothing to do
	}
        //- Read dictionary
        //void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "inertMultiComponentMixture.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
