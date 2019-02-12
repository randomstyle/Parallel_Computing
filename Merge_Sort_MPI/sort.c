/* assert */
#include <assert.h>
/* FILE, fopen, fclose, fscanf, rewind */
#include <stdio.h>
/* EXIT_SUCCESS, malloc, free, qsort */
#include <stdlib.h>
/* time, CLOCKS_PER_SEC */
#include <time.h>

static void
load(
  char const * const filename,
  int * const np,
  unsigned int ** const valsp
)
{
  int ret;
  int j, n;
  unsigned int dummy;
  FILE * fp;
  unsigned int * vals;

  /* open the file */
  fp = fopen(filename, "r");
  assert(fp);

  /* count the number of values in the file */
  for (n=0; fscanf(fp, "%u", &dummy)>0; ++n);
  assert(feof(fp));

  /* allocate memory for values */
  vals = malloc(n*sizeof(*vals));
  assert(vals);

  /* reset file pointer */
  rewind(fp);

  /* read in the values */
  for (j=0; j<n; ++j) {
    ret = fscanf(fp, "%u", &vals[j]);
    assert(1 == ret);
  }

  /* close file */
  ret = fclose(fp);
  assert(!ret);

  /* record output values */
  *np    = n;
  *valsp = vals;
}

static int
asc(
  void const * const a,
  void const * const b
)
{
  return (*(unsigned int*)a) < (*(unsigned int *)b) ? -1 : 1;
}

static void
sort_local(
  int const n,
  unsigned int * const vals
)
{
  qsort(vals, n, sizeof(*vals), asc);
}

static void
print_numbers(
  char const * const filename,
  int const n,
  unsigned int const * const vals
)
{
  int ret, i;
  FILE * fp;

  /* open file */
  fp = fopen(filename, "w");
  assert(fp);

  /* write list to fout */
  for (i=0; i<n; ++i) {
    fprintf(fp, "%u\n", vals[i]);
  }

  /* close file */
  ret = fclose(fp);
  assert(!ret);
}

static void
print_time(
  double const seconds
)
{
  printf("Sort Time: %0.04fs\n", seconds);
}

int
main(int argc, char ** argv)
{
  int n;
  clock_t s, e;
  unsigned int * vals;

  assert(argc > 2);

  load(argv[1], &n, &vals);

  s = clock(); /* use MPI_Wtime() in an MPI context */

  sort_local(n, vals);

  e = clock();

  print_numbers(argv[2], n, vals);
  print_time((double)(e-s)/CLOCKS_PER_SEC);

  free(vals);

  return EXIT_SUCCESS;
}
