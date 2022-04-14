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
    Foam::reactiveSurfaceConcentrationGradientFvPatchScalarField

Description
fixed gradient boundary conditions for explicit reaction rate

SourceFiles
    reactiveSurfaceConcentrationGradientFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef reactiveSurfaceConcentrationGradientFvPatchScalarFields_H
#define reactiveSurfaceConcentrationGradientFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class reactiveSurfaceConcentrationGradientFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class reactiveSurfaceConcentrationGradientFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data
    
    word omegaName_;


public:

    //- Runtime type information
    TypeName("reactiveSurfaceConcentrationGradient");


    // Constructors

        //- Construct from patch and internal field
        reactiveSurfaceConcentrationGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        reactiveSurfaceConcentrationGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  reactiveSurfaceConcentrationGradientFvPatchScalarField onto a new patch
        reactiveSurfaceConcentrationGradientFvPatchScalarField
        (
            const reactiveSurfaceConcentrationGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        reactiveSurfaceConcentrationGradientFvPatchScalarField
        (
            const reactiveSurfaceConcentrationGradientFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new reactiveSurfaceConcentrationGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        reactiveSurfaceConcentrationGradientFvPatchScalarField
        (
            const reactiveSurfaceConcentrationGradientFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new reactiveSurfaceConcentrationGradientFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
