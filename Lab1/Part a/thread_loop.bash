export OMP_NUM_THREADS=1
./test.ex

export OMP_NUM_THREADS=4
./test.ex

export OMP_NUM_THREADS=8
./test.ex

export OMP_NUM_THREADS=16
./test.ex

export OMP_NUM_THREADS=32
./test.ex

unset OMP_NUM_THREADS

