echo "Running hello.ex..."
mpirun -np 8 -hostfile hosts hello.ex
echo " "

echo "Running testB.ex..."
mpirun -np 8 -hostfile hosts testB.ex
echo " "

echo "Running testVec.ex..."
mpirun -np 8 -hostfile hosts testVec.ex
echo " "

echo "Running test_w.ex..."
mpirun -np 8 -hostfile hosts test_w.ex
echo " "

echo "Running sum_bcast.ex..."
mpirun -np 8 -hostfile hosts sum_bcast.ex
echo " "

echo "Running sum_bcast_inplace.ex..."
mpirun -np 8 -hostfile hosts sum_bcast_inplace.ex
echo " "

echo "Running sum_scatter.ex..."
mpirun -np 8 -hostfile hosts sum_scatter.ex
echo " "

echo "Running test_exp.ex..."
mpirun -np 8 -hostfile hosts test_exp.ex
echo " "

echo "Running test_exp_vec.ex..."
mpirun -np 8 -hostfile hosts test_exp_vec.ex
echo " "

echo "Running test_get_cnt.ex..."
mpirun -np 8 -hostfile hosts test_get_cnt.ex
echo " " 

echo "Running test_probe.ex..."
mpirun -np 8 -hostfile hosts test_probe.ex
echo " "

