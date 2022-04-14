
scalar meshCoNum = 0.0;
scalar meanMeshCoNum = 0.0;

if (mesh.nInternalFaces())
{
    surfaceScalarField SfUfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*mag(mesh.phi());

    meshCoNum = max(SfUfbyDelta/mesh.magSf())
        .value()*runTime.deltaT().value();

    meanMeshCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
        .value()*runTime.deltaT().value();
}

Info<< "Mesh Courant Number mean: " << meanMeshCoNum
    << " max: " << meshCoNum << endl;
