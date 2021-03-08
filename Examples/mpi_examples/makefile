all: hello.ex testB.ex testVec.ex test_w.ex

hello.ex: hello.c
	mpicc -o hello.ex hello.c

testB.ex: testB.c
	mpicc -o testB.ex testB.c

testVec.ex: testVec.c
	mpicc -o testVec.ex testVec.c

test_w.ex: test_w.c
	mpicc -o test_w.ex test_w.c

clean:
	rm *.ex


