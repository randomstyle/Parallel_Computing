#!/bin/bash
# 
#SBATCH --job-name=jtt
#SBATCH --output=jtt-srun.log
#SBATCH --nodes=5
#SBATCH --ntasks=140
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=1000

# export OMP_NUM_THREADS=10
export PATH=$PATH:/usr/lib/gcc/x86_64-redhat-linux/4.8.2/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/ohpc/pub/mpi/openmpi3-gnu7/3.0.0/lib

make floyd_mpi
mpirun -n 5 -N 1 ./floyd_mpi /home/011816337/213/project/Floyd-Parallel-Programming/matrices/1000.graph test.out
echo "==================check result==================="
diff test.out ori1000_test.out
echo "==================end==================="