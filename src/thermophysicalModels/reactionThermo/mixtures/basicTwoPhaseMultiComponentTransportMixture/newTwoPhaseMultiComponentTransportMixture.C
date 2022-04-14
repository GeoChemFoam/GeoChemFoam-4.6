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

#include "basicTwoPhaseMultiComponentTransportMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicTwoPhaseMultiComponentTransportMixture> Foam::basicTwoPhaseMultiComponentTransportMixture::New
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const volScalarField& alpha1
)
{
    word twoPhaseMultiComponentMixtureType;
    word MixtureTypeName1;
    word MixtureTypeName2;
    bool PhaseTransfer;

    // Enclose the creation of the thermophysicalProperties to ensure it is
    // deleted before the turbulenceModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermoPhysicalProperties",
                mesh.time().constant(),
                obj,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        
        PhaseTransfer= thermoDict.lookupOrDefault("PhaseTransfer",false);
        thermoDict.lookup("Phase1") >> MixtureTypeName1;
        thermoDict.lookup("Phase2") >> MixtureTypeName2;
    }

    Info<< "Selecting Mixture type " << MixtureTypeName1 << " for Phase 1" << endl;
    Info<< "Selecting Mixture type " << MixtureTypeName2 << " for Phase 2" << endl;

    if (PhaseTransfer==true)
    {
        twoPhaseMultiComponentMixtureType = "twoPhaseMultiComponentTransferMixture<" + MixtureTypeName1 + ',' + MixtureTypeName2 + '>';
    }
    else
    {
        twoPhaseMultiComponentMixtureType = "twoPhaseMultiComponentTransportMixture<" + MixtureTypeName1 + ',' + MixtureTypeName2 + '>';
    }
    
    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(twoPhaseMultiComponentMixtureType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("basicTwoPhaseMultiComponentTransportMixture::New(const fvMesh&, const objectRegistry&)")
            << "Unknown twoPhaseMultiComponentTransportMixture type "
            << twoPhaseMultiComponentMixtureType << nl << nl
            << "Valid twoPhaseMultiComponentTransportMixture types are:" << nl
            << fvMeshConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<basicTwoPhaseMultiComponentTransportMixture>(cstrIter()(mesh, obj,alpha1));
}

// ************************************************************************* //
