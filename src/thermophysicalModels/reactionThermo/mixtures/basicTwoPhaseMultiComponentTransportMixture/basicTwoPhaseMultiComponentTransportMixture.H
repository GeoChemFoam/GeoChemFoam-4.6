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
    Foam::basicTwoPhaseMultiComponentTransportMixture

Description
    Foam::basicTwoPhaseMultiComponentTransportMixture

SourceFiles
    basicTwoPhaseMultiComponentTransportMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicTwoPhaseMultiComponentTransportMixture_H
#define basicTwoPhaseMultiComponentTransportMixture_H

#include "basicMultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class basicTwoPhaseMultiComponentTransportMixture Declaration
\*---------------------------------------------------------------------------*/

class basicTwoPhaseMultiComponentTransportMixture
:
    public IOdictionary,
    public basicMultiComponentMixture
{

protected:

    // Protected data
        const fvMesh& mesh_;
        const volScalarField& alpha1_;

        string phiHScheme_;

        //- Diffusion and Henry's constant
        PtrList<dimensionedScalar> D1Y_;
        PtrList<dimensionedScalar> D2Y_;
        PtrList<dimensionedScalar> HY_;
        PtrList<dimensionedScalar> Mw_;

        //Total diffusion flux
        surfaceScalarField phiD_;

        //mass interfacial flux
        volScalarField Mflux_;

public:

    //- Runtime type information
    TypeName("basicTwoPhaseMultiComponentTransportMixture");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicTwoPhaseMultiComponentTransportMixture,
        fvMesh,
        (const fvMesh& mesh,
         const objectRegistry& obj,
         const volScalarField& alpha1
        ),
        (mesh, obj, alpha1)
    );


    // Constructors

        //- Construct from mesh and object registry
        basicTwoPhaseMultiComponentTransportMixture
        (
            const fvMesh&,
            const objectRegistry& obj,
            const volScalarField& alpha1
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<basicTwoPhaseMultiComponentTransportMixture> New(const fvMesh&, const objectRegistry&, const volScalarField& alpha1);

        static autoPtr<basicTwoPhaseMultiComponentTransportMixture> New(const fvMesh& mesh, const volScalarField& alpha1)
        {
            return New(mesh, mesh, alpha1);
        }


    //- Destructor
    virtual ~basicTwoPhaseMultiComponentTransportMixture();


    // Member functions

        // Access to phase parameters

        //- Return the species molecular diffusion list for phase 1
        const PtrList<dimensionedScalar>& D1Y() const
        {
            return D1Y_;
        }

        //- Return the species molecular diffusion for phase 1
        dimensionedScalar D1Y(const label i) const
        {
            return D1Y_[i];
        }

        //- Return the species molecular diffusion list for phase 2
        const PtrList<dimensionedScalar>& D2Y() const
        {
            return D2Y_;
        }

        //- Return the species molecular diffusion for phase 2
        dimensionedScalar D2Y(const label i) const
        {
            return D2Y_[i];
        }

        //- Return the species Henry's constant
        dimensionedScalar HY(const label i) const
        {
            return HY_[i];
        }

        //- Return molecular weight
        dimensionedScalar Mw(const label i) const
        {
            return Mw_[i];
        }

        //- Return mass interfacial flux
        const volScalarField Mflux() const
        {
            return Mflux_;
        }

        // calculate and return mean diffusion coefficient
        surfaceScalarField DmY(const label i) const
        {
            return fvc::interpolate(D1Y_[i]*alpha1_+HY_[i]*D2Y_[i]*(1-alpha1_))/fvc::interpolate(alpha1_+HY_[i]*(1-alpha1_));
        } 

        //calculate and return the henry transfer flux
        surfaceScalarField phiH(const label i) const;

        //calculate and return the upwind henry transfer flux
        surfaceScalarField phiHUp(const label i) const;

        //calculate and return the downwind henry transfer flux
        surfaceScalarField phiHDown(const label i) const;

        //calculate compression coefficient
        surfaceScalarField compressionCoeff(const label i);

        //- Update properties
        virtual void correct()=0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
