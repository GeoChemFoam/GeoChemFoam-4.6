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

Application
    multiSpeciesTransportFoam

Description
    Solves a transport equation with multiple passive species

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "multiComponentTransportMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    simpleControl simple(mesh);

#   include "createFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
			forAll(solutionSpecies, i)
			{
				volScalarField& Yi = speciesMixture.Y(i);
				dimensionedScalar DYi = speciesMixture.DY(i);

				fvScalarMatrix YiEqn
				(
					fvm::ddt(Yi)
					+ fvm::div(phi, Yi, "div(phi,Yi)")
					- fvm::laplacian(DYi, Yi)
				);

				YiEqn.solve(mesh.solutionDict().solver("Yi"));

		        Info<< solutionSpecies[i] << " concentration = "
			    << "  Min(Yi) = " << gMin(Yi)
			    << "  Max(Yi) = " << gMax(Yi)
			    << endl;
			}
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
