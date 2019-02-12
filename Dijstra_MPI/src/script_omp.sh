#!/bin/bash
# 
#SBATCH --job-name=pl2ap
#SBATCH --output=pl2ap-srun.log
# 
#SBATCH --ntasks=1
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=28
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=1000
# export OMP_NUM_THREADS=28
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# /home/david/programs/pl2ap/build/pl2ap pl2ap -t 0.9   
# /home/david/data/WikiWords200.csr
make dijkstra_omp
array=( 1 2 4 8 16 28 )
# array=( 20 )
for i in "${array[@]}"
do
  export OMP_NUM_THREADS=${i}
  echo ${i}
  /home/011816337/213/assignment3/src/dijkstra_omp /data/cmpe213sp18/pr3/data/100.graph 0 test.out ${i}
  /home/011816337/213/assignment3/src/dijkstra_omp /data/cmpe213sp18/pr3/data/500.graph 0 test.out ${i}
  /home/011816337/213/assignment3/src/dijkstra_omp /data/cmpe213sp18/pr3/data/1000.graph 0 test.out ${i}
  /home/011816337/213/assignment3/src/dijkstra_omp /data/cmpe213sp18/pr3/data/5000.graph 0 test.out ${i}
  /home/011816337/213/assignment3/src/dijkstra_omp /data/cmpe213sp18/pr3/data/10000.graph 0 test.out ${i}
  
done
echo "==================end==================="