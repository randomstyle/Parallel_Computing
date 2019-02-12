HEADERS := 

make:
	gcc-7 matmult2.c -fopenmp -o matmult2
	gcc matmult2.c -fopenmp -o matmult2
clean:
	-rm -f matmult2.o
	-rm -f matmult2

run:
	./matmult2 ./data/500-SQR-A.mat ./data/500-SQR-B.mat C.sol
	# ./matmult2 ./data/4000x1500-RCT-A.mat ./data/1500x3750-RCT-B.mat C.sol

run1:
	./matmult1 ./data/4000x1500-RCT-A.mat ./data/1500x3750-RCT-B.mat C2.sol


# ./matmult ./data/1000x2000-RCT-A.mat ./data/2000x5000-RCT-B.mat C.sol