CC=/sw/linux/mpi/gcc/openmpi/bin/mpic++
CC_RUN=/sw/linux/mpi/gcc/openmpi/bin/mpirun
target:
	${CC} -O3 --fast-math test_hist.cpp -o histo.out
	#${CC} -O3 --fast-math mpi_hist.cpp -o histo.out

run:
	${CC_RUN} -np 8 ./histo.out 8000000 20480
plot:
	python3 plot.py
clean:
	rm histo.out
