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
    Foam::solutionSurfaceMultiComponentMixture

Description
    Foam::solutionSurfaceMultiComponentMixture

SourceFiles
    solutionSurfaceMultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef solutionSurfaceMultiComponentMixture_H
#define solutionSurfaceMultiComponentMixture_H

#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solutionSurfaceMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

class solutionSurfaceMultiComponentMixture
:
    public basicMultiComponentMixture
{

protected:
    // Protected data

        //surface species mixture
        basicMultiComponentMixture surfaceMixture_;
        //surface species
        const speciesTable& surfaceSpecies_;
        //surface masters
        speciesTable surfaceMasters_;


    // Private member functions

        //- Construct as copy (not implemented)
        solutionSurfaceMultiComponentMixture(const solutionSurfaceMultiComponentMixture&);


public:


    // Constructors

        //- Construct from mesh and object registry
        solutionSurfaceMultiComponentMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&
        );

        //- Construct from mesh and object registry and phase name
        solutionSurfaceMultiComponentMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&,
            const word& phaseName
        );

        //- Construct from mesh
        solutionSurfaceMultiComponentMixture
        (
            const fvMesh&
        );


    // Member functions

        //- Return the masters species table
        basicMultiComponentMixture& surfaceMixture()
        {
            return surfaceMixture_;
        }

        const basicMultiComponentMixture& surfaceMixture() const
        {
            return surfaceMixture_;
        }

        //- Return the masters species table
        const speciesTable& surfaceMasters() const
        {
            return surfaceMasters_;
        }

        //- Read dictionary
        //void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "solutionSurfaceMultiComponentMixture.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
