#!/bin/bash

for type in normal cauchy
do
    for xopt in 0 1 2 10
    do
        for seed in 1024 4057 7090
        do
            sbatch run_sim.srun "-s ${seed}88${xopt} -d BURNIN=40000 -d NUMGENS=1000 -d N=10000 -d STEP=10 -d XOPT=${xopt}.0 ${type}.slim"
        done
    done
done
