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

#include "reactiveSurfaceConcentrationGradientFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactiveSurfaceConcentrationGradientFvPatchScalarField::
reactiveSurfaceConcentrationGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    omegaName_("omega")
{}


reactiveSurfaceConcentrationGradientFvPatchScalarField::
reactiveSurfaceConcentrationGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    omegaName_(dict.lookup("omega"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


reactiveSurfaceConcentrationGradientFvPatchScalarField::
reactiveSurfaceConcentrationGradientFvPatchScalarField
(
    const reactiveSurfaceConcentrationGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    omegaName_(ptf.omegaName_)
{}


reactiveSurfaceConcentrationGradientFvPatchScalarField::
reactiveSurfaceConcentrationGradientFvPatchScalarField
(
    const reactiveSurfaceConcentrationGradientFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    omegaName_(ptf.omegaName_)
{}


reactiveSurfaceConcentrationGradientFvPatchScalarField::
reactiveSurfaceConcentrationGradientFvPatchScalarField
(
    const reactiveSurfaceConcentrationGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    omegaName_(ptf.omegaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reactiveSurfaceConcentrationGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volScalarField& omega = this->db().objectRegistry::lookupObject<volScalarField>(omegaName_);
    const fvPatchField<scalar>& omegap = omega.boundaryField()[patch().index()];
    
    const volScalarField& R = this->db().objectRegistry::lookupObject<volScalarField>("R_"+omegaName_);
    const fvPatchField<scalar>& Rp = R.boundaryField()[patch().index()];
    
    const dictionary& speciesTransportProperties = db().lookupObject<IOdictionary> ("thermoPhysicalProperties");

    const dictionary& solutionSpeciesDict = speciesTransportProperties.subDict("solutionSpecies");
    const dictionary& kineticPhaseReactionDict = speciesTransportProperties.subDict("kineticPhaseReactions");
   
    const dictionary& subdict = solutionSpeciesDict.subDict(this->dimensionedInternalField().name());

    dimensionedScalar D(subdict.lookup("D"));
    
    const dictionary& kprSubDict = kineticPhaseReactionDict.subDict(omegaName_);
    const dictionary& kprSpeciesSubDict = kprSubDict.subDict("species").subDict(this->dimensionedInternalField().name());
    
    scalar scoeff(readScalar(kprSpeciesSubDict.lookup("scoeff")));
    
    gradient() = -scoeff*Rp*(1-omegap)/D.value();
  
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void reactiveSurfaceConcentrationGradientFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("omega") << omegaName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    reactiveSurfaceConcentrationGradientFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
