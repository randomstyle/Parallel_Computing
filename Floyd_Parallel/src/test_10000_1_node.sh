#!/bin/bash
#SBATCH --job-name=jessie
#SBATCH --output=jessie_1_node.log
#SBATCH --nodes=1
#SBATCH --ntasks=28


export PATH=$PATH:/usr/lib/gcc/x86_64-redhat-linux/4.8.2/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/ohpc/pub/mpi/openmpi3-gnu7/3.0.0/lib

make floyd_mpi_omp
echo "1 nodes"
array=( 5 10 20 25 )
for i in "${array[@]}"
do
echo ${i}
export OMP_NUM_THREADS=${i}
mpirun -N 1 ./floyd_mpi_omp /home/011816337/213/project/Floyd-Parallel-Programming/matrices/10000.mat /home/011816337/213/project/Floyd-Parallel-Programming/src/output/10000_mpi_1.out
done
echo "done"
