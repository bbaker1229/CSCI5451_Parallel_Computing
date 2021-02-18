module load soft/cuda
nvcc testSaxpy.cu
echo "Running with only computations:"
./a.out
nvcc testSaxpyC.cu
echo "Running with computations and data transfer:"
./a.out
rm a.out
