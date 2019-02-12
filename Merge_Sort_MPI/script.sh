#!/bin/bash
# 
#SBATCH --job-name=pl2ap
#SBATCH --output=pl2ap-srun.log
# 
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=4000
# export OMP_NUM_THREADS=28
# export OMP_PLACES=cores
# export OMP_PROC_BIND=spread

make sort_mpi
# nodes=( 1 2 5 10 15 20 )
# procs=( 1 2 4 8 16 28 )
# procs=( 28 )
# array=( 28 )

# for j in "${procs[@]}"
# do
    # a=$(( 1*j ))
    # # echo "$a"
# mpirun /home/011816337/213/assignment2/sort_mpi /data/cmpe213sp18/pr2/data/sm.dat result.dat
# mpirun -np 28 /home/011816337/213/assignment2/sort_mpi /data/cmpe213sp18/pr2/data/md.dat result.dat
    # mpirun /home/011816337/213/assignment2/sort_mpi /data/cmpe213sp18/pr2/data/lg.dat result.dat
mpirun -np 10 /home/011816337/213/assignment2/sort_mpi /data/cmpe213sp18/pr2/data/lg.dat result.dat
# done


echo "==================end==================="