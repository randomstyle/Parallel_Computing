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
#include "omp.h"
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
    a = malloc(n*n*sizeof(*a));
    assert(a);

    int size;
    int rem = n%num_proc;

    if (my_rank >= num_proc-rem && rem != 0) size = n/num_proc + 1;
    else size = n/num_proc;
    
    a = calloc(n*size, sizeof(*a));
    assert(a);

    if (my_rank == 0){
        int size1;
        for (int k=0; k<num_proc; ++k){
            if (k >= num_proc-rem && rem != 0) size1 = n/num_proc + 1;
            else size1 = n/num_proc;
            
            a_temp = calloc(n*size1, sizeof(*a_temp));
            assert(a_temp);

            for (i=0; i<size1; ++i) {
                for (j=0; j<n; ++j) {
                    if (k == 0){
                        ret = fscanf(fp, "%f", &a(i,j));
                        assert(1 == ret);
                        if(a(i,j) == 0.0){
                            a(i,j) = INF;
                        }
                    }
                    else{
                        ret = fscanf(fp, "%f", &a_temp(i,j));
                        assert(1 == ret);
                        if(a_temp(i,j) == 0.0){
                            a_temp(i,j) = INF;
                        }
                    }
                }
            }
        
            if (k != 0) MPI_Send(a_temp, n*size1, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
            free(a_temp);
        }
    }
    else MPI_Recv(a, n*size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

void floydWarshall (int const n, float * a, float ** res_p, int my_rank, int num_proc, int size) {
    int i, j, k;
    float * res;
    float * row_k;
    double tts, tte;
    row_k = malloc(n * sizeof(*row_k));
    
    if(my_rank == 0){
        res = malloc(n*n*sizeof(*res));
    }

    // if(my_rank == 1){
    //     for (int t1 = 0; t1 < size; t1++) {
	// 	    for (int t2 = 0; t2 < n; t2++){
	// 		    if (a(t1, t2) ==0.0) {
	// 			    printf("t1: %d, t2 %d, a %f \n", t1, t2, a(t1, t2));
	// 		    }
    //         }
	//     }
    // }

    // omp_set_num_threads(20);
    // #pragma omp parallel
    // {
    //     printf("abcd %d\n", omp_get_num_threads());
    // }

    for (k = 0; k < n; k++) {
        int root = k/size;
        if (my_rank == root) {
            int k_temp = k-size*my_rank;
			for (int p = 0; p < n; p++) {
				row_k[p] = a(k_temp,p);
			}
		}
		MPI_Bcast(row_k, n, MPI_FLOAT, root, MPI_COMM_WORLD);
        if (my_rank == 0) tts = MPI_Wtime();
        #pragma omp parallel for
        for (i = 0; i < size; i++) {
            float a_i_k = a(i,k);
            for (j = 0; j < n; j++) {
                float temp = a_i_k + row_k[j];
                if (temp < a(i,j))
                    a(i,j) = temp;
            }
        }
        if (my_rank == 0) {tte = MPI_Wtime();printf("inner check time: %f", tte-tts);}
    }
    MPI_Gather(a, size*n, MPI_FLOAT, res, size*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    *res_p = res;
}

int main(int argc, char ** argv) {
    int n;
    float * a;
    float * res;
    double ts, te, ts_o, te_o;
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
    /* 
        make a cartesian (block of processes communicator) this is for 2D partition
    */

    // MPI_Comm cart_com;
    // int ndims = 2;
    // int reorder = FALSE;
    // int dim_size[2] = {3, 4};
    // int periods[2] = {TRUE, TRUE};

    // int coord[2];

    // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    // MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &cart_com);

    // if(my_rank == 2){
    //     MPI_Cart_coords(cart_com, my_rank, 2, coord);
    //     printf("%d %d\n", coord[0], coord[1]);
    // }
    // // printf("This is the %d\n", my_rank);

    /* 
        read data
    */
    if (my_rank ==0){ 
        printf("new world size %d\n", num_proc);
        #pragma omp parallel
        {
            if(omp_get_thread_num() == 0) printf("thread num %d\n", omp_get_num_threads());
        }
    }
    load(argv[1], &n, &a, my_rank, num_proc, &size);

    MPI_Barrier(MPI_COMM_WORLD);
    ts = MPI_Wtime();
    ts_o = omp_get_wtime();
    floydWarshall(n, a, &res, my_rank, num_proc, size);
    te = MPI_Wtime();
    te_o = omp_get_wtime();

    //printf("My rank is: %d Operation Time: %0.04fs\n", my_rank, te-ts);
    if(my_rank == 0){
        printf("Operation Time OMP: %0.04fs\n", te_o-ts_o);
        printf("Operation Time MPI: %0.04fs\n", te-ts);
        print_numbers(argv[2],n,res);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}