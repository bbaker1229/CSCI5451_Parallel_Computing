echo "Running hello.ex..."
mpirun -np 8 -hostfile hosts hello.ex

echo "Running testB.ex..."
mpirun -np 8 -hostfile hosts testB.ex

echo "Running testVec.ex..."
mpirun -np 8 -hostfile hosts testVec.ex

echo "Running test_w.ex..."
mpirun -np 8 -hostfile hosts test_w.ex

