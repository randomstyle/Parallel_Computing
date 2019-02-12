/* errno */
#include <errno.h>

/* fopen, fscanf, fprintf, fclose */
#include <stdio.h>

/* EXIT_SUCCESS, EXIT_FAILURE, malloc, free */
#include <stdlib.h>
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


static int mult_mat(size_t const n, size_t const m, size_t const p,
                    double const * const A, double const * const B,
                    double ** const Cp)
{
  size_t i, j, k;
  double sum;
  double * C=NULL;

  if (!(C=malloc(n*p*sizeof(*C)))) {
    goto cleanup;
  }

  for (i=0; i<n; ++i) {
    for (j=0; j<p; ++j) {
      for (k=0,sum=0.0; k<m; ++k) {
        sum += A[i*m+k] * B[k*p+j];
      }
      C[i*p+j] = sum;
    }
  }

  *Cp = C;

  return 0;
  
  cleanup:
  free(C);

  /*failure:*/
  return -1;
}


int main(int argc, char * argv[])
{
  // size_t stored an unsigned integer
  size_t n, m, p, mm;
  double * A=NULL, * B=NULL, * C=NULL, start, end;

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
  if (mult_mat(n, m, p, A, B, &C)) {
    perror("error");
    goto failure;
  }
  end = omp_get_wtime();
  fprintf(stderr, "Time take: %f", end-start);

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
