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

#include "collapsingFace.H"
#include "polyMesh.H"
#include "faceSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(collapsingFace, 0);

addToRunTimeSelectionTable(topoSetSource, collapsingFace, word);

addToRunTimeSelectionTable(topoSetSource, collapsingFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::collapsingFace::usage_
(
    collapsingFace::typeName,
    "\n    Usage: collapsingFace <tol>\n\n"
    "    Select faces which collapse within tol to internalFaces\n"
);

// Construct from Istream
bool Foam::collapsingFace::checkHit(label faceI, label faceJ) const
{
	bool hit = false;

	const pointField& meshPoints = mesh_.points();

	{
		//Info << "check hit " << faceI << " " << faceJ << endl;
		point faceCentre = mesh_.faces()[faceI].centre(meshPoints);

		vector n = mesh_.faceAreas()[faceI];
		vector v = -n;


		pointHit nearest = mesh_.faces()[faceJ].ray(faceCentre, v, meshPoints, intersection::FULL_RAY);

		if (nearest.distance() < 0 || nearest.distance() > tol_) nearest.setMiss(true);
		hit = nearest.hit();
	}

	if (!hit)
	{
		//Info << "check hit " << faceI << " " << faceJ << endl;
		point faceCentre = mesh_.faces()[faceJ].centre(meshPoints);

		vector n = mesh_.faceAreas()[faceJ];
		vector v = -n;


		pointHit nearest = mesh_.faces()[faceI].ray(faceCentre, v, meshPoints, intersection::FULL_RAY);
		if (nearest.distance() < 0 || nearest.distance() > tol_) nearest.setMiss(true);
		hit = nearest.hit();
	}
	return hit;
}

// Construct from Istream
bool Foam::collapsingFace::checkIntersect(label faceI, label faceJ) const
{
	bool intersect=false;

	const pointField& meshPoints = mesh_.points();

	{
		const pointField& facePoint = mesh_.faces()[faceJ].points(meshPoints);
		point faceCentre = mesh_.faces()[faceJ].centre(meshPoints);

		forAll(facePoint, pointi)
		{
			vector n = facePoint[pointi] - faceCentre;
			scalar length = mag(n);
			vector v = n / length;
			pointHit nearestFace = mesh_.faces()[faceI].ray(faceCentre, v, meshPoints, intersection::FULL_RAY);
			if (nearestFace.hit() && nearestFace.distance() > 0 && nearestFace.distance() < 0.999*length)
			{
				intersect = true;
				break;
			}
			forAll(facePoint, pointj)
			{	
				if (pointj == pointi)
				{
					continue;
				}
				n = facePoint[pointi] - facePoint[pointj];
				length = mag(n);
				v = n / length;
				nearestFace = mesh_.faces()[faceI].ray(facePoint[pointj], v, meshPoints, intersection::FULL_RAY);
				if (nearestFace.hit() && nearestFace.distance() > 0 && nearestFace.distance() < 0.999*length)
				{
					intersect = true;
					break;
				}
			}
			if (intersect) break;
		}
	}

	return intersect;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::collapsingFace::collapsingFace
(
    const polyMesh& mesh,
	const scalar tol
)
:
    topoSetSource(mesh),
	tol_(tol)
{
}


// Construct from dictionary
Foam::collapsingFace::collapsingFace(const polyMesh& mesh, const dictionary& dict)
:
    topoSetSource(mesh),
	tol_(readScalar(dict.lookup("tol")))
{
}


// Construct from Istream
Foam::collapsingFace::collapsingFace(const polyMesh& mesh, Istream& is)
:
    topoSetSource(mesh),
	tol_(readScalar(checkIs(is)))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::collapsingFace::~collapsingFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::collapsingFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Collapsing face with new or add not implemented " << endl;

    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing faces which collapse "
		<< " (within " << 0 << " and " << tol_ << " to internal faces) ..."  << endl;

        dynamicLabelList toBeRemoved(set.size()/10);
        dynamicLabelList intersectingFaces(set.size()/10);

        forAllIter(topoSet, set, iter1)
        {
            label faceI = iter1.key();
			bool remove =false;

			forAllIter(topoSet, set, iter2)
			{
                           
			    label faceJ = iter2.key();
				if (faceJ==faceI)
				{
					continue;
				}

			 	scalar d = mag(mesh_.faceCentres()[faceI] - mesh_.faceCentres()[faceJ]);
				if (d > 10 * tol_)
				{
					continue;
				}

				if (!checkHit(faceI, faceJ))
				{
					continue;
				}

				if (checkIntersect(faceI, faceJ))
				{
					break;
				}

				remove = true;
				{
					break;
				}
			}


			if (remove)
			{
				forAllIter(topoSet, set, iter2)
				{
					label faceJ = iter2.key();
					if (faceJ==faceI)
					{
						continue;
					}

					scalar d = mag(mesh_.faceCentres()[faceI] - mesh_.faceCentres()[faceJ]);
					if (d > 10 * tol_)
					{
						continue;
					}

					if (checkIntersect(faceI, faceJ))
					{
						remove = false;
						break;
					}
				}
			}
							

			if (remove) toBeRemoved.append(faceI);
        }

        forAll(toBeRemoved, i)
        {
            set.erase(toBeRemoved[i]);
        }
    }
}


// ************************************************************************* //
