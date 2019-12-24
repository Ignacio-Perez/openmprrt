#!/bin/bash

OUT=cputimes.csv

REPETITIONS=30

TIMEFORMAT='%3R'

BIN=rrt

MAPS="1 2 3 4"

ALG="1 2"

CPUS="1 2 4"

echo CPUs ";" ALG ";" MAP ";" TIME > $OUT

for c in $CPUS
do
	export OMP_NUM_THREADS=$c
	for a in $ALG
	do
		for m in $MAPS
		do
			for r in `seq 1 $REPETITIONS`
			do
				t=`(time ./$BIN $m $a 1) 2>&1 >/dev/null`
				echo $c $a $m $t
				echo $c ";" $a ";" $m ";" $t >> $OUT
			done
		done
	done
done

rm temp
