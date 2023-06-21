#!/bin/bash

for type in normal cauchy
do
    for xopt in 0.0 1.0 2.0 10.0
    do
        for seed in 123 456 789
        do
            sbatch run_sim.srun "-s ${seed} -d BURNIN=40000 -d NUMGENS=10000 -d N=10000 -d XOPT=${xopt} ${type}.slim"
        done
    done
done
