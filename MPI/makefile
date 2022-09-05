
CC=mpic++

target:
	${CC} -O3 --fast-math histogram.cpp -o histo.out
run:
	mpirun -np 2 ./histo.out 8000000 20480
clean:
	rm histo.out
