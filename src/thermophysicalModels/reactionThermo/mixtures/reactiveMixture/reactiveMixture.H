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
    Foam::reactiveMixture

Description
    Foam::reactiveMixture

SourceFiles
    reactiveMixture.C

\*---------------------------------------------------------------------------*/

#ifndef reactiveMixture_H
#define reactiveMixture_H

#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class reactiveMixture Declaration
\*---------------------------------------------------------------------------*/

class reactiveMixture
:
    public basicMultiComponentMixture
{

protected:
    // Private member functions

        wordList kineticPhase_;
        wordList kineticPhaseReaction_;
        
        PtrList<dimensionedScalar> Ke_;
        PtrList<dimensionedScalar> k0_;
        PtrList<scalarList> scoeff_;
        PtrList<labelList> kiSpeciesIndex_;
        PtrList<PtrList<dimensionedScalar> > ki_;
        
        PtrList<volScalarField> Omega_;
        PtrList<volScalarField> R_;
        
        //index mapping to index of 1st si component
	int* component1Index_;
	//index mapping to index of 2nd si component
	int* component2Index_;
        
        
        //- Construct as copy (not implemented)
        reactiveMixture(const reactiveMixture&);


public:


    // Constructors

        //- Construct from mesh and object registry
        reactiveMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&
        );

        //- Construct from mesh and object registry and phase name
        reactiveMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const objectRegistry&,
            const word&
        );

        //- Construct from mesh
        reactiveMixture
        (
            const fvMesh&
        );


    // Destructor
	virtual ~reactiveMixture();
	
    // Member functions
    void correct();

    void setSaturation(const volScalarField& alpha)
	{
        //nothing to do
	}
	
    const wordList& kineticPhases()
    {
        return kineticPhase_;
    }

    const wordList& kineticPhaseReactions()
    {
        return kineticPhaseReaction_;
    }
    
    const volScalarField& omega(label i)
    {
        return Omega_[i];
    }
    
    const scalarList& scoeff(label i)
    {
        return scoeff_[i];
    }
    
    const volScalarField& Rate(label i)
    {
        return R_[i];
    }       
        //- Read dictionary
        //void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "reactiveMixture.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
