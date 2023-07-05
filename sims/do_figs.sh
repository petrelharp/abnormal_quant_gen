#!/bin/bash
# OLD VERSION: takes a loooong time to run the selected sims
# see make_paper_figs.sh

for type in normal cauchy
do
    for seed in 8135 8791
    do
        sbatch run_sim.srun "-s ${seed} -d BURNIN=40000 -d NUMGENS=10000 -d N=10000 neutral_${type}.slim"
        for xopt in 0 1 2 10
        do
            sbatch run_sim.srun "-s ${seed}55${xopt} -d BURNIN=40000 -d NUMGENS=10000 -d N=10000 -d STEP=10 -d XOPT=${xopt}.0 ${type}.slim"
        done
    done
done
