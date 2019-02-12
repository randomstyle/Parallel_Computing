SRC := sort_mpi.c

default: sort_mpi

sort: sort.c
	gcc -O3 -Wall -Wextra -o $@ sort.c

sort_mpi: $(SRC)
	mpicc -O3 -Wall -Wextra -o $@ $(SRC)

clean:
	-rm -f sort sort_mpi 
