#!/bin/bash

echo "SLOTS PROCS MEAN SIGMA" > results.txt

# Hostfile2
for i in 2 4 8
do 
  echo -n "2 ${i} " >> results.txt; mpirun -np ${i} -hostfile hostfile2 main_Stats.ex | awk '{print $2 " " $4}' >> results.txt
done

# Hostfile4
for i in 2 4 8 16
do
  echo -n "4 ${i} " >> results.txt; mpirun -np ${i} -hostfile hostfile4 main_Stats.ex | awk '{print $2 " " $4}' >> results.txt
done

# Hostfile8
for i in 2 4 8 16 32
do 
  echo -n "8 ${i} " >> results.txt; mpirun -np ${i} -hostfile hostfile8 main_Stats.ex | awk '{print $2 " " $4}' >> results.txt
done

# Hostfile16
for i in 2 4 8 16 32
do 
  echo -n "16 ${i} " >> results.txt; mpirun -np ${i} -hostfile hostfile16 main_Stats.ex | awk '{print $2 " " $4}' >> results.txt
done

# Hostfile32
for i in 2 4 8 16 32
do
  echo -n "32 ${i} " >> results.txt; mpirun -np ${i} -hostfile hostfile32 main_Stats.ex | awk '{print $2 " " $4}' >> results.txt
done
