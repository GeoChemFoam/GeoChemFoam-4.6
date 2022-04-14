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
    Foam::twoPhaseMultiComponentTransferMixture

Description
    Foam::twoPhaseMultiComponentTransferMixture

SourceFiles
    twoPhaseMultiComponentTransferMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMultiComponentTransferMixture_H
#define twoPhaseMultiComponentTransferMixture_H

#include "basicTwoPhaseMultiComponentTransportMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class twoPhaseMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template <class MixtureType1,class MixtureType2>
class twoPhaseMultiComponentTransferMixture
:
    public basicTwoPhaseMultiComponentTransportMixture
{
    // Private data
        static const word phase1Name;
        static const word phase2Name;

        //- Species data
        MixtureType1 phase1SpeciesMixture_;
        MixtureType2 phase2SpeciesMixture_;


    // Private member functions

        //- Construct as copy (not implemented)
        twoPhaseMultiComponentTransferMixture(const twoPhaseMultiComponentTransferMixture&);


public:

    //- Runtime type information
    TypeName("twoPhaseMultiComponentTransferMixture");

    // Constructors

        //- Construct from mesh
        twoPhaseMultiComponentTransferMixture
        (
            const fvMesh&,
            const objectRegistry&,
            const volScalarField& alpha1
        );

        //- Construct from mesh
        twoPhaseMultiComponentTransferMixture
        (
            const fvMesh&,
            const volScalarField& alpha1
        );

    // Member functions

        //correct interfacial transfer
        void correct();

        //- Read dictionary
        //void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "twoPhaseMultiComponentTransferMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
