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
#include <omp.h>

#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
/* define access directions for matrices */
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]
#define a_res(R,C) a_res[ROWMJR(R,C,ln,n)]
 
#define INF 99999.0
 
void floydWarshall (int const n, float * a, float ** const res_p) {

    int i, j, k;

    for (k = 0; k < n; k++) {
        #pragma omp parallel for
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                float temp = a(i,k) + a(k,j);
                if (temp < a(i,j))
                    a(i,j) = temp;
            }
        }
    }
    *res_p = a;
}

static void load(char const * const filename, int * const np, float ** const ap) {
    int i, j, n, ret;
    FILE * fp=NULL;
    float * a;

    /* open the file */
    fp = fopen(filename, "r");
    assert(fp);

    /* get the number of nodes in the graph */
    ret = fscanf(fp, "%d", &n);
    assert(1 == ret);

    /* allocate memory for local values */
    a = malloc(n*n*sizeof(*a));
    assert(a);

    /* read in roots local values */
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j) {
            ret = fscanf(fp, "%f", &a(i,j));
            assert(1 == ret);
            if(a(i,j) == 0.0){
                a(i,j) = INF;    
            }
        }
    }

    /* close file */
    ret = fclose(fp);
    assert(!ret);

    /* record output values */
    *np = n;
    *ap = a;
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

// driver program to test above function
int main(int argc, char ** argv) {   
    int n;
    float * a;
    float * res;
    double ts, te;

    if(argc < 3){
        printf("Invalid number of arguments.\nUsage: dijkstra <graph> <output_file>.\n");
        return EXIT_FAILURE;
    }

    load(argv[1], &n, &a);

    ts = omp_get_wtime();
    floydWarshall(n, a, &res);
    te = omp_get_wtime();

    print_numbers(argv[2],n,res);
    printf("Operation Time: %0.04fs\n", te-ts);

    return 0;
}