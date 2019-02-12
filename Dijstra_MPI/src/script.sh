#!/bin/bash
# 
#SBATCH --job-name=pl2ap
#SBATCH --output=pl2ap-srun.log
# 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


#SBATCH --time=10:00
#SBATCH --mem-per-cpu=4000

make dijkstra_mpi


mpirun -np 1 /home/011816337/213/assignment3/src/dijkstra_mpi /data/cmpe213sp18/pr3/data/100.graph 0 test.out 1 1

# 1,2,5,10  1,2,4,8,16,28
# 1,2,4,8,16,28
# 2,4,8,16,32,64
# 5,10,20,40,80,140
# 10,20,40,80,160,280

# 1,2,4,5,8,10,16,20,28,32,40,64,80,140,160,280
# 100     1,2,4,5,10,20,25,50
# 500     1,2,4,5,10,20,25,50,100,125
# 1000    1,2,4,5,10,20,25,40,50,100,125,200,250
# 5000    1,2,4,5,10,20,25,40,50,100,125,200,250
# 10000   1,2,4,5,10,20,25,40,50,80,100,125,200,250
echo "==================end==================="
