#!/bin/bash
for ((i=0;i<=6;i++))
do
    date >> finish.deck
    mv Dataa$[i*20+20]qe Data
    ./bin/gztrace
    mv Data Dataa$[i*20+20]qe
    echo Dataa$[i*20+20]qe is finished! >> finish.deck
#    echo Dataa$[ia*2+250] is running! >> finish.deck
#    echo Dataa$[ia*2+250] | /public/software/mpi/openmpi/1.8.5/intel/bin/mpirun -np $NP -machinefile $PBS_NODEFILE ./bin/epoch2d 
done
