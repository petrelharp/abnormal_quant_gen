#!/bin/bash
SLIM="sbatch run_sim.srun"
# SLIM=slim

$SLIM -s 831 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 -d 'EPSILON=c(0,0.2,1.0,4.0,10)' neutral_cauchy.slim
$SLIM -s 831 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 -d 'EPSILON=c(0,0.2,1.0,4.0,10)' neutral_normal.slim

$SLIM -s 707 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 -d STEP=20 -d 'EPSILON=c(0,0.2,1.0,4.0,10)' cauchy.slim
$SLIM -s 707 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 -d STEP=20 -d 'EPSILON=c(0,0.2,1.0,4.0,10)' normal.slim

$SLIM -s 8135 -d BURNIN=40000 -d NUMGENS=1000 -d N=10000 -d 'EPSILON=c(0,0.2,1.0,4.0,10)' neutral_normal.slim
$SLIM -s 8135 -d BURNIN=40000 -d NUMGENS=1000 -d N=10000 -d 'EPSILON=c(0,0.2,1.0,4.0,10)' neutral_cauchy.slim

$SLIM -s 8136 -d BURNIN=40000 -d NUMGENS=1000 -d N=10000 -d 'EPSILON=c(0,0.2,1.0,4.0,10)' neutral_stable.slim
