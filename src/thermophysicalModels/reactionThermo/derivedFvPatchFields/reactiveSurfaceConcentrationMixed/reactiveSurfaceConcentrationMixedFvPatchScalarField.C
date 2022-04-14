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

#include "reactiveSurfaceConcentrationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactiveSurfaceConcentrationMixedFvPatchScalarField::reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    k_(0.0),
    scoeff_(0.0)
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


reactiveSurfaceConcentrationMixedFvPatchScalarField::reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const reactiveSurfaceConcentrationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    k_(ptf.k_),
    scoeff_(ptf.scoeff_)
{}


reactiveSurfaceConcentrationMixedFvPatchScalarField::reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    k_(readScalar(dict.lookup("k"))),
    scoeff_(readScalar(dict.lookup("scoeff")))
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}


reactiveSurfaceConcentrationMixedFvPatchScalarField::reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const reactiveSurfaceConcentrationMixedFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    k_(tppsf.k_),
    scoeff_(tppsf.scoeff_)
{}


reactiveSurfaceConcentrationMixedFvPatchScalarField::reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const reactiveSurfaceConcentrationMixedFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    k_(tppsf.k_),
    scoeff_(tppsf.scoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reactiveSurfaceConcentrationMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const dictionary& speciesTransportProperties = db().lookupObject<IOdictionary> ("thermoPhysicalProperties");

    const dictionary& solutionSpeciesDict = speciesTransportProperties.subDict("solutionSpecies");

    const dictionary& subdict = solutionSpeciesDict.subDict(this->dimensionedInternalField().name());

    dimensionedScalar D(subdict.lookup("D"));


    scalarField lambda = scoeff_*k_/D.value()/(this->patch().deltaCoeffs());

    valueFraction() = lambda/(lambda + 1);


    mixedFvPatchScalarField::updateCoeffs();
}

void reactiveSurfaceConcentrationMixedFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
    os.writeKeyword("scoeff") << scoeff_ << token::END_STATEMENT << nl;
    mixedFvPatchScalarField::write(os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, reactiveSurfaceConcentrationMixedFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
