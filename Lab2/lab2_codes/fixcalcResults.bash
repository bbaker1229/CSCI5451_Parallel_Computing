#!/bin/bash

echo "SLOTS PROCS MEAN SIGMA" > fix_results.txt

# Hostfile2
for i in 2 4 8
do 
  echo -n "2 ${i} " >> fix_results.txt; mpirun -np ${i} -hostfile hostfile2 fix_main_Stats.ex | awk '{print $2 " " $4}' >> fix_results.txt
done

# Hostfile4
for i in 2 4 8 16
do
  echo -n "4 ${i} " >> fix_results.txt; mpirun -np ${i} -hostfile hostfile4 fix_main_Stats.ex | awk '{print $2 " " $4}' >> fix_results.txt
done

# Hostfile8
for i in 2 4 8 16 32
do 
  echo -n "8 ${i} " >> fix_results.txt; mpirun -np ${i} -hostfile hostfile8 fix_main_Stats.ex | awk '{print $2 " " $4}' >> fix_results.txt
done

# Hostfile16
for i in 2 4 8 16 32
do 
  echo -n "16 ${i} " >> fix_results.txt; mpirun -np ${i} -hostfile hostfile16 fix_main_Stats.ex | awk '{print $2 " " $4}' >> fix_results.txt
done

# Hostfile32
for i in 2 4 8 16 32
do
  echo -n "32 ${i} " >> fix_results.txt; mpirun -np ${i} -hostfile hostfile32 fix_main_Stats.ex | awk '{print $2 " " $4}' >> fix_results.txt
done

echo "NODES PPR MEAN SIGMA" > fix_ppr_results.txt

for i in 1 2 3 4
do
for j in 2 4 8 16
do
  let "val = ${i} * ${j}"
  echo -n "${i} ${j} " >> fix_ppr_results.txt; mpirun -np $val -hostfile hostfile64 -map-by ppr:${j}:node fix_main_Stats.ex | awk '{print $2 " " $4}' >> fix_ppr_results.txt
done
done

