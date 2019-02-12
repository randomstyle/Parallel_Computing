/* errno */
#include <errno.h>

/* fopen, fscanf, fprintf, fclose */
#include <stdio.h>

/* EXIT_SUCCESS, EXIT_FAILURE, malloc, free */
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

static int load_mat(char const * const fname, size_t * const np,
                    size_t * const mp, double ** const matp)
{
  size_t i, j, n, m;
  double * mat=NULL;
  FILE * fp;

  if (!(fp=fopen(fname, "r"))) {
    goto failure;
  }

  if (2 != fscanf(fp, "%zu %zu", &n, &m)) {
    goto cleanup;
  }

  if (!(mat=malloc(n*m*sizeof(*mat)))) {
    goto cleanup;
  }

  for (i=0; i<n; ++i) {
    for (j=0; j<m; ++j) {
      if (feof(fp)) {
        goto cleanup;
      }
      fscanf(fp, "%lf", mat+i*m+j);
    }
  }

  if (fclose(fp)) {
    goto failure;
  }

  *np   = n;
  *mp   = m;
  *matp = mat;

  return 0;

  cleanup:
  free(mat);
  fclose(fp);

  failure:
  return -1;
}

static int save_mat(char const * const fname, size_t const n,
                    size_t const m, double const * const mat)
{
  size_t i, j;
  FILE * fp;

  if (!(fp=fopen(fname, "w"))) {
    goto failure;
  }

  fprintf(fp, "%zu %zu\n", n, m);

  for (i=0; i<n; ++i) {
    for (j=0; j<m; ++j) {
      fprintf(fp, "%10.4lf ", mat[i*m+j]);
    }
    fprintf(fp, "\n");
  }

  if (fclose(fp)) {
    goto failure;
  }

  return 0;

  failure:
  return -1;
}


static int mult_mat(double const * const A, double const * const B, size_t const n, size_t const p, size_t const m, size_t const rstart, size_t const rend, size_t const cstart, size_t const cend,size_t const qstart, size_t const qend, double *C)
{
  
  size_t i, j, k;
  double sum;

  for (i = rstart; i <= rend; i++) {
    for (j = cstart; j <= cend; j++) {
      // if (qstart == 0) 
      sum = 0.0;
      for (k = qstart; k <= qend; k++) {
        // fprintf(stderr, "%f\n",A[i*m+k]);
        // fprintf(stderr, "%f\n",B[k*p+j]);
        // sleep(1);
        sum += A[i*m+k] * B[k*p+j];
      } /* for q */

      C[i*p+j] += sum;
    } /* for c */
  } /* for r */
  
  return 0;
  
}

static int matrix_matrix_mult_by_tiling (double** dst, double* src1, double* src2, size_t nr, size_t nc, size_t nq, size_t rtilesize, size_t ctilesize, size_t qtilesize){ 
  /* matrix_matrix_mult_by_tiling */
  
  size_t rstart, rend, cstart, cend, qstart, qend;
  double * C=NULL;
  if (!(C=malloc(nr*nc*sizeof(*C)))) {
    goto cleanup;
  }
  #pragma omp parallel
  {
  #pragma omp for private(rstart, rend, cstart, cend, qstart, qend)
  for (rstart = 0; rstart < nr; rstart += rtilesize) {
    rend = rstart + rtilesize - 1;
    if (rend >= nr) rend = nr - 1;
    for (cstart = 0; cstart < nc; cstart += ctilesize) {
      cend = cstart + ctilesize - 1;
      if (cend >= nc) cend = nc - 1;
      for (qstart = 0; qstart < nq; qstart += qtilesize) {
        qend = qstart + qtilesize - 1;
        if (qend >= nq) qend = nq - 1;
        mult_mat(src1, src2, nr, nc, nq, rstart, rend, cstart, cend, qstart, qend, C);
        
      } /* for qstart */
    } /* for cstart */
  } /* for rstart */
  }
  *dst = C;
  
  return 0;

  cleanup:
  free(C);
  return -1;
} /* matrix_matrix_mult_by_tiling */

int main(int argc, char * argv[])
{
  // size_t stored an unsigned integer
  size_t n, m, p, mm, rtilesize, ctilesize, qtilesize;
  rtilesize = 10;
  ctilesize = 10;
  qtilesize = 10;
  double * A=NULL, * B=NULL, * C=NULL, start, end;;

  if (argc != 4) {
    fprintf(stderr, "usage: matmult A.mat B.mat C.sol\n");
    goto failure;
  }

  if (load_mat(argv[1], &n, &m, &A)) {
    perror("error");
    goto failure;
  }

  if (load_mat(argv[2], &mm, &p, &B)) {
    perror("error");
    goto failure;
  }

  if (m != mm) {
    fprintf(stderr, "dimensions do not match: %zu x %zu, %zu x %zu\n",
      n, m, mm, p);
  }
  start = omp_get_wtime();
  if (matrix_matrix_mult_by_tiling(&C, A, B, n, p, m, rtilesize, ctilesize, qtilesize)) {
    perror("error");
    goto failure;
  }
  end = omp_get_wtime();
  fprintf(stderr, "Time take: %f\n", end-start);

  if (save_mat(argv[3], n, p, C)) {
    perror("error");
    goto failure;
  }

  free(A);
  free(B);
  free(C);

  return EXIT_SUCCESS;

  failure:
  free(A);
  free(B);
  free(C);
  return EXIT_FAILURE;
}
