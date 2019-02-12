FLOYD := floyd_mpi_omp.c
FLOYD_SERIAL := floyd.c
FLOYD_OMP := floyd_omp.c
FLOYD_MPI := floyd_mpi.c


default: floyd

floyd_mpi_omp: $(FLOYD)
	mpicc -fopenmp -O3 -Wall -Wextra -o $@ $<

floyd: $(FLOYD_SERIAL)
	gcc -O3 -Wall -Wextra -o $@ $<

floyd_omp: $(FLOYD_OMP)
	gcc -fopenmp -O3 -Wall -Wextra -o $@ $<

floyd_mpi: $(FLOYD_MPI)
	mpicc -O3 -Wall -Wextra -o $@ $<

clean: 
	rm -f floyd floyd_omp floyd_mpi_omp a.out floyd_mpi
