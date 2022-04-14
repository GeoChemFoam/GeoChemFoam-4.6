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

#include "reactingWallFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
	surfaceMasters_(dict.subDict("surfaceMasters").toc()),
	density_(surfaceMasters_.size())
{	
	forAll(surfaceMasters_, i)
	{
		const dictionary& surfaceMastersDict = dict.subDict("surfaceMasters");
		const dictionary& subdict = surfaceMastersDict.subDict(surfaceMasters_[i]);
		density_[i] = readScalar(subdict.lookup("density"));
	}
}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const reactingWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
	surfaceMasters_(ptf.surfaceMasters_),
	density_(ptf.density_)
{
}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const reactingWallFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
	surfaceMasters_(ptf.surfaceMasters_),
	density_(ptf.density_)
{}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const reactingWallFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
	surfaceMasters_(ptf.surfaceMasters_),
	density_(ptf.density_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reactingWallFvPatchScalarField::write(Ostream& os) const
{
	dictionary dict;

	forAll(surfaceMasters_, i)
	{
		dictionary subdict;
		subdict.add("density",density_[i]);
		dict.add(surfaceMasters_[i],subdict);
	}

	os.writeKeyword("surfaceMasters") << dict << endl;
    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    reactingWallFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
