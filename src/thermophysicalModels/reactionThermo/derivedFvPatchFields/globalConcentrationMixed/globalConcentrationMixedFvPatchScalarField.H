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
    Foam::globalConcentrationMixedFvPatchScalarField

Description
 Mixed boundary conditions for concentration gradient on a wetted wall

SeeAlso
    Foam::mixedFvPatchScalarField

SourceFiles
    globalConcentrationMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef globalConcentrationMixedFvPatchScalarField_H
#define globalConcentrationMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class globalConcentrationMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class globalConcentrationMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private data

    // Henry's coefficient
    scalar H_;

public:

    //- Runtime type information
    TypeName("globalConcentrationMixed");


    // Constructors

        //- Construct from patch and internal field
        globalConcentrationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        globalConcentrationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  globalConcentrationMixedFvPatchScalarField onto a new patch
        globalConcentrationMixedFvPatchScalarField
        (
            const globalConcentrationMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        globalConcentrationMixedFvPatchScalarField
        (
            const globalConcentrationMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new globalConcentrationMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        globalConcentrationMixedFvPatchScalarField
        (
            const globalConcentrationMixedFvPatchScalarField&,
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
                new globalConcentrationMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions



        // Evaluation functions

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
