if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxMeshCo/(meshCoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaT().value(),
            maxDeltaT
        )
    );

    scalar newTime = runTime.value()+runTime.deltaT().value();
    if (newTime>runTime.endTime().value()) runTime.setDeltaT(runTime.endTime().value()-runTime.value(),false); 
    else if (newTime+0.5*runTime.deltaT().value() > runTime.endTime().value()) runTime.setDeltaT((runTime.endTime().value()-runTime.value())/2.0,false); 
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //
