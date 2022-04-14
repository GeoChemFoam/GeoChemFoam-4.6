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

InClass
    Foam::twoPhaseMultiComponentMixture

Description

\*---------------------------------------------------------------------------*/

#ifndef makeTwoPhaseMultiComponentMixture_H
#define makeTwoPhaseMultiComponentMixture_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeTwoPhaseMultiComponentTransportMixture(TwoPhaseType,MixtureType1,MixtureType2)                                      \
                                                                              \
typedef TwoPhaseType<MixtureType1,MixtureType2>                               \
    TwoPhaseType##MixtureType1##MixtureType2;                                 \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    TwoPhaseType##MixtureType1##MixtureType2,                                 \
    #TwoPhaseType"<"#MixtureType1","#MixtureType2">",                         \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicTwoPhaseMultiComponentTransportMixture,                              \
    TwoPhaseType##MixtureType1##MixtureType2,                                 \
    fvMesh                                                                    \
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
