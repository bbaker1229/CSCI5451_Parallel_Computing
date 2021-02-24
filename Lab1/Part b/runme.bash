module load soft/cuda
nvcc testSaxpy.cu
echo "Running with only computations:" > testSaxpy.sav
./a.out >> testSaxpy.sav
nvcc testSaxpyC.cu
echo "Running with computations and data transfer:" > testSaxpyC.sav
./a.out >> testSaxpyC.sav
rm a.out
