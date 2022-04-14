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
    Foam::interfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{

private:
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

		//- smoothing coefficient for alphaS
		scalar cSK_;
		//-number of smoothing cycle for alphaS
		label nSK_;
		//- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
		scalar cPc_;
		//-number of non-orthogonal corrector loop
		label nNonOrthogonalCorrectors_;
		/*----------------------------------------------*/

        //- Surface tension
        dimensionedScalar sigma_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField nHatf_;
        volScalarField K_;

		//interface normal vector at cell center
		volVectorField nI_;
		//interface curvature at face center
        surfaceScalarField Kf_;
		//capillary pressure
		volScalarField pc_;
		label pcRefCell_;
		scalar pcRefValue_;
		//Capillary flow	
		surfaceScalarField phic_;

		//capillary pressure residual
		scalar eqnResidual_;

		/*----------------------------------------------*/

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfaceProperties(const interfaceProperties&);
        void operator=(const interfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::GeometricBoundaryField& nHat,
            surfaceVectorField::GeometricBoundaryField& gradAlphaf
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();

		//Re-calculate capillary flux
		void calculatePhic();

public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        const volScalarField& K() const
        {
            return K_;
        }

        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

        tmp<volScalarField> sigmaK() const
        {
            return sigma_*K_;
        }

		//function correct now include calculate capillary flux
        void correct()
        {
            calculateK();
            calculatePhic();
        }

		//get capillary pressure
		const volScalarField& pc()
		{
			return pc_;
		}


		//get capillary flow of pc
		const surfaceScalarField& phic()
		{
			return phic_;
		}
		/*----------------------------------------------*/


		//get capillary pressure equation residual
		const scalar& eqnResidual()
		{
			return eqnResidual_;
		}

                volScalarField nearInterface() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
