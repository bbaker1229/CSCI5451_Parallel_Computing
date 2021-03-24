all: hello.ex testB.ex testVec.ex test_w.ex sum_bcast.ex sum_bcast_inplace.ex sum_scatter.ex test_exp.ex test_exp_vec.ex test_get_cnt.ex test_probe.ex

hello.ex: hello.c
	mpicc -o hello.ex hello.c

testB.ex: testB.c
	mpicc -o testB.ex testB.c

testVec.ex: testVec.c
	mpicc -o testVec.ex testVec.c

test_w.ex: test_w.c
	mpicc -o test_w.ex test_w.c

sum_bcast.ex: sum_bcast.c
	mpicc -o sum_bcast.ex sum_bcast.c

sum_bcast_inplace.ex: sum_bcast_IN_Place.c
	mpicc -o sum_bcast_inplace.ex sum_bcast_IN_Place.c

sum_scatter.ex: sum_scatter.c
	mpicc -o sum_scatter.ex sum_scatter.c

test_exp.ex: test_exp.c
	mpicc -o test_exp.ex test_exp.c

test_exp_vec.ex: test_exp_vec.c
	mpicc -o test_exp_vec.ex test_exp_vec.c

test_get_cnt.ex: test_get_cnt.c
	mpicc -o test_get_cnt.ex test_get_cnt.c -lm

test_probe.ex: test_probe.c
	mpicc -o test_probe.ex test_probe.c -lm

clean:
	rm *.ex


