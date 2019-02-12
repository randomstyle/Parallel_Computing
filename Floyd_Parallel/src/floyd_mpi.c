/* assert */
#include <assert.h>
/* INFINITY */
#include <math.h>
/* FILE, fopen, fclose, fscanf, rewind */
#include <stdio.h>
/* EXIT_SUCCESS, malloc, calloc, free */
#include <stdlib.h>
/* time, CLOCKS_PER_SEC */
#include <time.h>
#include <mpi.h>
// #include <omp.h>

#define TRUE 1
#define FALSE 0
#define INF 99999.0

#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]
#define a_res(R,C) a_res[ROWMJR(R,C,ln,n)]
#define a_temp(R,C) a_temp[ROWMJR(R,C,ln,n)]
#define d(R,C) d[ROWMJR(R,C,ln,n)]


static void load(char const * const filename, int * const np, float ** const ap, int my_rank, int num_proc, int * const size_p) {
    int i, j, n, ret;
    FILE * fp=NULL;
    float * a;
    float * a_temp;

    /* open the file */
    fp = fopen(filename, "r");
    assert(fp);

    /* get the number of nodes in the graph */
    ret = fscanf(fp, "%d", &n);
    assert(1 == ret);

    /* allocate memory for local values */
    // a = malloc(n*n*sizeof(*a));
    // assert(a);

    //arfaerfearf
    int size;
    int rem = n%num_proc;

    if (my_rank >= num_proc-rem && rem != 0) size = n/num_proc + 1;
    else size = n/num_proc;
    
    a = malloc(n*size*sizeof(*a));
    assert(a);

    if (my_rank == 0){
        int size1;
        for (int k=0; k<num_proc; ++k){
            if (k >= num_proc-rem && rem != 0) size1 = n/num_proc + 1;
            else size1 = n/num_proc;
            // printf("%d", size1);
            a_temp = malloc(n*size1*sizeof(*a_temp));
            assert(a_temp);

            for (i=0; i<size1; ++i) {
                for (j=0; j<n; ++j) {
                    if (k == 0){
                        ret = fscanf(fp, "%f", &a(i,j));
                        assert(1 == ret);
                        if(a(i,j) == 0){
                            a(i,j) = INFINITY;
                        }
                        // printf("a is %f\n", a(i,j));
                    } 
                    else{
                        ret = fscanf(fp, "%f", &a_temp(i,j));
                        assert(1 == ret);
                        if(a_temp(i,j) == 0){
                            a_temp(i,j) = INFINITY;
                        }
                        // printf("a_temp is %f\n", a_temp(i,j));
                    }
                }
            }
        
            if (k != 0) MPI_Send(a_temp, n*size1, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
            free(a_temp);
        }
    }
    else MPI_Recv(a, n*size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   
    /* read in roots local values */
    // for (i=0; i<n; ++i) {
    //     for (j=0; j<n; ++j) {

    //     }
    // }

    /* close file */
    ret = fclose(fp);
    assert(!ret);

    /* record output values */
    *np = n;
    *ap = a;
    *size_p = size;
}



static void print_numbers(char const * const filename, int const n, float const * const a) {
    int i, j;
    FILE * fout;

    /* open file */
    if(NULL == (fout = fopen(filename, "w"))) {
        fprintf(stderr, "error opening '%s'\n", filename);
        abort();
    }

    /* write numbers to fout */
    for(i=0; i<n; ++i) {
        for(j=0; j<n; ++j){
            if(i==j) fprintf(fout, "%10.4s\t", "inf");
            else fprintf(fout, "%10.4f\t", a(i,j));
        }
        fprintf(fout, "%s", "\n");
    }

    fclose(fout);
}

void floydWarshall (int const n, float * a, float ** res_p, int my_rank, int size) {

    //size is number of rows for each process

    int i, j, k;
	float * row_k;
	row_k = malloc(n*sizeof(*row_k));

    for (k = 0; k < n; k++) {
        // Pick all vertices as source one by one
       
		//one to all broad cast 	
		int root = k / size;

	//	if (my_rank == 0) printf("k:%d,  k_new: %d, root: %d\n", k, k_new, root);
		if (my_rank == root) {
			
			int k_temp = k - size * my_rank;
		
			for (int p = 0; p < n; p++) {
				
				row_k[p] = a(k_temp, p);

				
			}
		}
		
		
		MPI_Bcast(row_k, n, MPI_FLOAT, root, MPI_COMM_WORLD);
	
		
	
        for (i = 0; i < size; i++) {
            // Pick all vertices as destination for the above picked source
			float a_i_k = a(i, k);
            for (j = 0; j < n; j++) {
                // If vertex k is on the shortest path from i to j, then update the value of dist[i][j]
				float temp = a_i_k + row_k[j];
                //float temp = a(i,k) + a(k,j);
				if (temp < a(i, j)) {
					a(i, j) = temp;					
				}
            }
        }
    }
		
	float * global_a;
	global_a = malloc(n * n * sizeof(*global_a));
	MPI_Gather(a, size*n, MPI_FLOAT, global_a, size*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
    *res_p = global_a;
}

int main(int argc, char ** argv) {
    int n;
    float * a;
    float * res;
    double ts, te;
    int my_rank;
    int num_proc;
    int size;

    if(argc < 3){
        printf("Invalid number of arguments.\nUsage: dijkstra <graph> <output_file>.\n");
        return EXIT_FAILURE;
    }

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    
    load(argv[1], &n, &a, my_rank, num_proc, &size);
	    
	MPI_Barrier(MPI_COMM_WORLD);
	ts = MPI_Wtime();
    floydWarshall(n, a, &res, my_rank, size);
    te = MPI_Wtime();

	//printf("rank %d, Operation Time: %0.04fs\n", my_rank, te - ts);
	if (my_rank == 0) {
		printf("Number of processes: %d, Operation Time: %0.04fs\n", num_proc, te - ts);
		print_numbers(argv[2], n, res);
	}    

    MPI_Finalize();
    return EXIT_SUCCESS;
}