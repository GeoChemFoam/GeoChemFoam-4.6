#decompose
echo -e "decomposePar"
decomposePar > decomposePar.out

# run simpleDBSFoam
echo -e "Run simpleDBSFoam in parallel"
mpiexec -np 4 simpleDBSFoam -parallel > simpleDBSFoam.out
echo -e "reconstructPar"
reconstructPar -latestTime > reconstructPar.out

rm -rf processor*

echo -e "processPoroPerm" 
processPoroPerm > processPoroPerm.out

rm *.out
