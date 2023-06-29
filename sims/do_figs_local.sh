#!/bin/bash

if [ $# -lt 1 ]
then
    echo "Usage: $0 seed"
fi

seed=$1

for xopt in 0 1 2 10
do
    for type in normal cauchy
    do
        echo "starting at $(date)" >> slurm.log
        bash run_sim.srun "-s ${seed}99${xopt} -d BURNIN=10000 -d NUMGENS=1000 -d N=1000 -d XOPT=${xopt}.0 ${type}.slim"
        echo "ending: $(date)" >> slurm.log
    done
done
