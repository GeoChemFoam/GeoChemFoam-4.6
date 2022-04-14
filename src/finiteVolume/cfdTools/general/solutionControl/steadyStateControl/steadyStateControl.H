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
    Foam::steadyStateControl

Description
    SIMPLE control class to supply convergence information/checks for
    the SIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef steadyStateControl_H
#define steadyStateControl_H

#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class steadyStateControl Declaration
\*---------------------------------------------------------------------------*/

class steadyStateControl
:
    public solutionControl
{

protected:

    // Protected Data

        //- Initialised flag
        bool initialised_;

       label iter_counter;

       label nCorr_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        void read();

        //- Disallow default bitwise copy construct
        steadyStateControl(const steadyStateControl&);

        //- Disallow default bitwise assignment
        void operator=(const steadyStateControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("steadyStateControl");


    // Constructors

        //- Construct from mesh
        steadyStateControl(fvMesh& mesh);


    //- Destructor
    virtual ~steadyStateControl();


    // Member Functions

        // Solution control
//- Return true if all convergence checks are satisfied
        bool criteriaSatisfied();

            //- Loop loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
