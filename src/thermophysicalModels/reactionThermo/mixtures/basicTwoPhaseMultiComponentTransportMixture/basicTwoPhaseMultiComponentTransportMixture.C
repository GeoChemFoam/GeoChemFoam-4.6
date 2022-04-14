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
#include "upwind.H"
#include "downwind.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicTwoPhaseMultiComponentTransportMixture, 0);
    defineRunTimeSelectionTable(basicTwoPhaseMultiComponentTransportMixture, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicTwoPhaseMultiComponentTransportMixture::basicTwoPhaseMultiComponentTransportMixture
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const volScalarField& alpha1
)
:
    IOdictionary
    (
        IOobject
        (
            "thermoPhysicalProperties",
            mesh.time().constant(),
            obj,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    basicMultiComponentMixture(*this, this->subDict("solutionSpecies").toc(), mesh, obj),
    mesh_(mesh),
    alpha1_(alpha1),
	D1Y_(this->subDict("solutionSpecies").toc().size()),
	D2Y_(this->subDict("solutionSpecies").toc().size()),
	HY_(this->subDict("solutionSpecies").toc().size()),
    Mw_(this->subDict("solutionSpecies").toc().size()),
    phiD_
	(
		IOobject
		(
			"phiD",
            alpha1_.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("phiD", dimMass/dimTime, 0.0)
	),
    Mflux_
	(
		IOobject
		(
			"Mflux",
            alpha1_.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("MFlux", dimMass/dimVolume/dimTime, 0.0)
	)
{

    //******************************************************************//
    // scheme for phiH
    //******************************************************************//

	OStringStream buf;

	const tokenList& tokens = mesh_.schemesDict().divScheme("div(phiH,Yi)");

	label len = tokens.size();

	forAll (tokens, tok)
	{
        buf << mesh_.schemesDict().divScheme("div(phiH,Yi)")[tok];

		if (--len)
		{
	        buf << ' ';
	    }
	}
  
    phiHScheme_ = buf.str();

	forAll(species_, i)
	{
	    const dictionary& solutionSpeciesDict = this->subDict("solutionSpecies");
		const dictionary& subdict = solutionSpeciesDict.subDict(species_[i]);
		D1Y_.set
		(
			i,
			new dimensionedScalar(subdict.lookup("D1"))
		);
		D2Y_.set
		(
			i,
			new dimensionedScalar(subdict.lookup("D2"))
		);
		HY_.set
		(
			i,
			new dimensionedScalar(subdict.lookup("H"))
		);
		Mw_.set
		(
			i,
			new dimensionedScalar(subdict.lookup("Mw"))
		);
	}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicTwoPhaseMultiComponentTransportMixture::~basicTwoPhaseMultiComponentTransportMixture()
{}

// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //
//calculate and return the henry transfer flux
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture::phiH(const label i) const
{
    return DmY(i)*(1-HY_[i])/fvc::interpolate((alpha1_+HY_[i]*(1-alpha1_)))*fvc::snGrad(alpha1_)*mesh_.magSf();
}

//calculate and return the upwind henry transfer flux
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture::phiHUp(const label i) const
{
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)*mesh_.magSf();
	surfaceScalarField alphaUp = upwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
    return DmY(i)/(alphaUp+(1-alphaUp)*HY_[i])*fvc::snGrad(alpha1_)*mesh_.magSf();
}

//calculate and return the downwind henry transfer flux
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture:: phiHDown(const label i) const
{
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)*mesh_.magSf();
	surfaceScalarField alphaDown = downwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
    return -HY_[i]*DmY(i)/(alphaDown+(1-alphaDown)*HY_[i])*fvc::snGrad(alpha1_)*mesh_.magSf();

}

//calculate compression coeeficient for species i
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture::compressionCoeff(const label i)
{
    //Direction of interfacial flux
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)*mesh_.magSf();

    //upwind and downwind alpha1
	surfaceScalarField alphaUp = upwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
	surfaceScalarField alphaDown = downwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);

    //upwind and downwnd Yi
    surfaceScalarField YiUp   = upwind<scalar>(mesh_,fluxDir).interpolate(Y_[i]);
    surfaceScalarField YiDown = downwind<scalar>(mesh_,fluxDir).interpolate(Y_[i]);
    
    dimensionedScalar sgn = sign(max(alphaDown*YiDown)-max((1-alphaUp)*YiUp));

    //normal compression coefficient
    surfaceScalarField deltaYi1 = max(-max(Y_[i]),min(max(Y_[i]),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
    
    //standard compression coefficient
    surfaceScalarField deltaYi2 = max(-max(Y_[i]),min(max(Y_[i]),YiDown/(alphaDown+(1-alphaDown)*HY_[i])-HY_[i]*YiUp/(alphaUp+(1-alphaUp)*HY_[i])));

    return sgn*max(mag(deltaYi1),mag(deltaYi2));

}
