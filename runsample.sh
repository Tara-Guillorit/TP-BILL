#!/bin/bash

for i in $(seq 1 9);
do
	echo "running $i"
	srun python3 src/sample_stats.py -f 0.05 -d 150 -ec 10000 270000 -lc 2000 -s 1 -sa $i -it "1-5" -o results/sample_stats/$i
done	
