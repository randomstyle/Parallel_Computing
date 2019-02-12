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
export OMP_NUM_THREADS=28
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# /home/david/programs/pl2ap/build/pl2ap pl2ap -t 0.9   
# /home/david/data/WikiWords200.csr
gcc matmult2.c -fopenmp -o matmult
# array=( 1 4 8 12 16 20 24 28 )
array=( 28 )
for i in "${array[@]}"
do
  export OMP_NUM_THREADS=${i}
  echo $OMP_NUM_THREADS ": "
  /home/011816337/213/assignment1/matmult /data/cmpe213sp18/pr1/data/1000-SQR-A.mat /data/cmpe213sp18/pr1/data/1000-SQR-B.mat C.sol
  /home/011816337/213/assignment1/matmult /data/cmpe213sp18/pr1/data/2000-SQR-A.mat /data/cmpe213sp18/pr1/data/2000-SQR-B.mat C.sol
  /home/011816337/213/assignment1/matmult /data/cmpe213sp18/pr1/data/1000x2000-RCT-A.mat /data/cmpe213sp18/pr1/data/2000x5000-RCT-B.mat C.sol
  /home/011816337/213/assignment1/matmult /data/cmpe213sp18/pr1/data/4000x1500-RCT-A.mat /data/cmpe213sp18/pr1/data/1500x3750-RCT-B.mat C.sol
done
echo "==================end==================="