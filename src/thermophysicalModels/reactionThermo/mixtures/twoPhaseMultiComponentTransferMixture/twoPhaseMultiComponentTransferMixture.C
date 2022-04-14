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

#include "twoPhaseMultiComponentTransferMixture.H"
#include "upwind.H"
#include "downwind.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class MixtureType1, class MixtureType2>
const word Foam::twoPhaseMultiComponentTransferMixture<MixtureType1,MixtureType2>::phase1Name = "phase1";

template<class MixtureType1, class MixtureType2>
const word Foam::twoPhaseMultiComponentTransferMixture<MixtureType1,MixtureType2>::phase2Name = "phase2";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType1, class MixtureType2>
Foam::twoPhaseMultiComponentTransferMixture<MixtureType1,MixtureType2>::twoPhaseMultiComponentTransferMixture
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const volScalarField& alpha1
)
:
    basicTwoPhaseMultiComponentTransportMixture(mesh, obj, alpha1),
	phase1SpeciesMixture_(*this, this->subDict("solutionSpecies").toc(), mesh, obj, phase1Name),
	phase2SpeciesMixture_(*this, this->subDict("solutionSpecies").toc(), mesh, obj, phase2Name)
{
}

template<class MixtureType1, class MixtureType2>
Foam::twoPhaseMultiComponentTransferMixture<MixtureType1,MixtureType2>::twoPhaseMultiComponentTransferMixture
(
    const fvMesh& mesh,
    const volScalarField& alpha1
)
:
    basicTwoPhaseMultiComponentTransportMixture(mesh, mesh, alpha1),
	phase1SpeciesMixture_(*this, this->subDict("solutionSpecies").toc(), mesh, mesh, phase1Name),
	phase2SpeciesMixture_(*this, this->subDict("solutionSpecies").toc(), mesh, mesh, phase2Name)
{
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
template<class MixtureType1, class MixtureType2>
void Foam::twoPhaseMultiComponentTransferMixture<MixtureType1,MixtureType2>::correct()
{
    //calculate alpha downwind
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)*mesh_.magSf();
    surfaceScalarField alphaUp = upwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
    surfaceScalarField alphaDown = downwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
    surfaceScalarField alphaf    = fvc::interpolate(alpha1_);

    //sum of phase mass fraction
    surfaceScalarField sumX = 1 - max(alphaf, alphaDown-alphaUp);

    //Re-initialize transfer flux
    phiD_ *=0;

    if (phiHScheme_ == "Gauss upwind")
    {
		forAll(species_, i)
		{
			phiD_+=Mw_[i]/(1-sumX+1e-4)*
				 (
				 	DmY(i)*fvc::snGrad(Y_[i])*mesh_.magSf()
		           -fvc::flux(phiHUp(i),Y_[i],"div(phiH,Yi)")
		           -fvc::flux(phiHDown(i),Y_[i],"div(phiH,Yi)")
				 );
		}
    }
    else if (phiHScheme_ == "Gauss linear")
    {
		forAll(species_, i)
		{
			phiD_+=Mw_[i]/(1-sumX+1e-4)*
				 (
				 	DmY(i)*fvc::snGrad(Y_[i])*mesh_.magSf()
				   -fvc::flux(phiH(i),Y_[i],"div(phiH,Yi)")
				 );
		}
    }
    else
    {
        Info<< "div(phiH,Yi) should be equal to Gauss linear or Gauss upwind"
		<< endl
		<< abort(FatalError);
    }

    //Compute flux
	Mflux_ = fvc::div(phiD_*alphaDown)-alpha1_*fvc::div(phiD_);

    volScalarField alpha2 = 1-alpha1_;

	//compute Yi1 and Yi2
	forAll(species_, i)
	{
		volScalarField& Yi = Y_[i];
		volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
		volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
		dimensionedScalar HYi = HY_[i];
		Y1i = Yi/(alpha1_+HYi*(1-alpha1_));
		Y2i = HYi*Yi/(alpha1_+HYi*(1-alpha1_));
	}

    //set saturation
    phase1SpeciesMixture_.setSaturation(alpha1_);
    phase2SpeciesMixture_.setSaturation(alpha2);

    //correct each phase
    phase1SpeciesMixture_.correct();
    phase2SpeciesMixture_.correct();

	//compute Yi from Y1i and Y2i
	forAll(species_, i)
	{
		volScalarField& Yi = Y_[i];
		volScalarField& Y1i = phase1SpeciesMixture_.Y(i);
		volScalarField& Y2i = phase2SpeciesMixture_.Y(i);
		Yi = Y1i*alpha1_ + Y2i*(1-alpha1_);
	}
    
}
