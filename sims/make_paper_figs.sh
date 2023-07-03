slim -s 831 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 neutral_cauchy.slim
slim -s 831 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 neutral_normal.slim

slim -s 707 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 -d STEP=20 cauchy.slim
slim -s 707 -d BURNIN=20000 -d NUMGENS=1000 -d N=1000 -d STEP=20 normal.slim
