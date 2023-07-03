#!/bin/bash

for type in normal cauchy
do
    for xopt in 0 1 2 10
    do
        for seed in 123 456 789
        do
            sbatch run_sim.srun "-s ${seed}${xopt} -d BURNIN=40000 -d NUMGENS=10000 -d N=10000 -d XOPT=${xopt}.0 ${type}.slim"
        done
    done
done