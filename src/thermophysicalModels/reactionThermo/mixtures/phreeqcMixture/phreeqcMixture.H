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
    Foam::phreeqcMixture

Description
    Foam::phreeqcMixture

SourceFiles
    phreeqcMixture.C

\*---------------------------------------------------------------------------*/

#ifndef phreeqcMixture_H
#define phreeqcMixture_H

#include "solutionSurfaceMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PhreeqcMixture Declaration
\*---------------------------------------------------------------------------*/

class phreeqcMixture
:
    public solutionSurfaceMultiComponentMixture
{
    // Private data
	//initialise phreeqc module 
	void initialise();

	//OpenFOAM mesh
	const fvMesh& mesh_;

	//surface 
	volScalarField Surf_;

    //ionic strength
    volScalarField I_;

    //surface potential
    volScalarField psi_;

	//index of phreeqc module
	int id_;

	//water saturation
	double* saturation_;

	//index mapping from solution species to phreeqc component
	int* componentSolutionIndex_;
	//index mapping from surface species to phreeqc component
	int* componentSurfaceIndex_;
	//concentration
	double* concentration_;
	//surface concentration
	double* surfConcentration_;
	//surface area
	double* surfArea_;
	//surface potential
	double* surfPotential_;

	//selected output
	double* selectedOutput_;


    // Private member functions

        //- Construct as copy (not implemented)
        phreeqcMixture(const phreeqcMixture&);


public:


    // Constructors

        //- Construct from mesh and object registry
        phreeqcMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&
        );

        //- Construct from mesh and object registry and phase name
        phreeqcMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&,
            const word&
        );

        //- Construct from mesh
        phreeqcMixture
        (
            const fvMesh&
        );


    // Destructor
	virtual ~phreeqcMixture();

    // Member functions

	//solve reaction step function
	void correct();

    void setSaturation(const volScalarField& alpha);

        //- Read dictionary
        //void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "phreeqcMixture.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
