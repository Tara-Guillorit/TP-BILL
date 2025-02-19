#!/bin/bash

for i in $(seq 1 5)
do
	echo "parsing $i"
	srun python3 src/iteration_stats.py -f 0.15 -d 400 -ec 20000 270000 -lc 3000 -s 1 -sa 1-10 -it $i -o results/iteration_stats/v3.0/P$i
done
