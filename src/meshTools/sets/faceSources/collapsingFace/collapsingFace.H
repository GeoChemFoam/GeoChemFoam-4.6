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
    Foam::collapsingFace

Description
    A topoSetSource to select faces that collapse into internal faces.

SourceFiles
    collapsingFace.C

\*---------------------------------------------------------------------------*/

#ifndef collapsingFace_H
#define collapsingFace_H

#include "topoSetSource.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class collapsingFace Declaration
\*---------------------------------------------------------------------------*/

class collapsingFace
:
    public topoSetSource
{

private:

        //- Add usage string
        static addToUsageTable usage_;

        //- Tolerance (i.e. cos of angle between normal_ and faceNormal)
        const scalar tol_;

		bool checkHit(label faceI, label faceJ) const;

		bool checkIntersect(label faceI,label faceJ) const;



public:

    //- Runtime type information
    TypeName("collapsingFace");

    // Constructors

        //- Construct from components
        collapsingFace
        (
            const polyMesh& mesh,
            const scalar tol
        );

        //- Construct from dictionary
        collapsingFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        collapsingFace(const polyMesh& mesh, Istream&);


    // Destructor

        virtual ~collapsingFace();


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
