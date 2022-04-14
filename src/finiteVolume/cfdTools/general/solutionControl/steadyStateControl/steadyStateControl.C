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

#include "steadyStateControl.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyStateControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::steadyStateControl::read()
{
    solutionControl::read(true);
    const dictionary& steadyStateDict =dict();
    nCorr_ = steadyStateDict.lookupOrDefault<label>("nCorrectors",1000);
}


bool Foam::steadyStateControl::criteriaSatisfied()
{
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.solutionDict().solverPerformanceDict();
    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldI = applyToField(variableName);

        if (fieldI != -1)
        {
            scalar lastResidual = 0;
            const scalar residual =
                maxResidual(variableName, iter().stream(), lastResidual);

            checked = true;

            bool absCheck = lastResidual < residualControl_[fieldI].absTol;
            achieved = achieved && absCheck;

            if (debug)
            {
                Info<< algorithmName_ << " solution statistics:" << endl;

                Info<< "    " << variableName << ": tolerance = " << residual
                    << " (" << residualControl_[fieldI].absTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyStateControl::steadyStateControl(fvMesh& mesh)
:
    solutionControl(mesh, "STEADYSTATE"),
    initialised_(false),
    iter_counter(0)
{
    read();

    Info<< nl;

    if (residualControl_.empty())
    {
        Info<< algorithmName_ << ": no convergence criteria found. "
            << "Calculations will run for " << mesh_.time().endTime().value()
            << " steps." << nl << endl;
    }
    else
    {
        Info<< algorithmName_ << ": convergence criteria" << nl;
        forAll(residualControl_, i)
        {
            Info<< "    field " << residualControl_[i].name << token::TAB
                << " tolerance " << residualControl_[i].absTol
                << nl;
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::steadyStateControl::~steadyStateControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::steadyStateControl::loop()
{
    read();

    if (initialised_)
    {
	iter_counter++;
        if (iter_counter==nCorr_ + 1)
        {
            Info << nl << algorithmName_ << ": not converged within "
                 << nCorr_ << " iterations" << endl;
            return false;
        }
        if (criteriaSatisfied())
        {
            Info<< nl << algorithmName_ << " solution converged in "
                << iter_counter << " iterations" << nl << endl;

            // Set to finalise calculation
            return false;
        }
        else
        {
            storePrevIterFields();
        }
    }
    else
    {
        initialised_ = true;
        storePrevIterFields();
    }


    return true;
}


// ************************************************************************* //
