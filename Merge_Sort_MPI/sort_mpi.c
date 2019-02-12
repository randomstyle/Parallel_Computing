/* assert */
#include <assert.h>
/* FILE, fopen, fclose, fscanf, rewind */
#include <stdio.h>
/* EXIT_SUCCESS, malloc, free, qsort */
#include <stdlib.h>
/* time, CLOCKS_PER_SEC */
#include <mpi.h>

static int ryan_binary_search(int l, int r, unsigned int * values, unsigned int val){
  int ori_l = l;
  while (l <= r){
    int m = l+(r-l)/2;

    if (values[m] == val)
        return m-ori_l;

    if (values[m] < val){
        l = m + 1;
    }
    else{
        r = m - 1;
    }
  }

  return l - ori_l;
}
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
  int n, rem;
  double s, e;
  unsigned int * my_data;
  unsigned int * vals;
  int * my_size;
  int * displ;
  int * my_datas;

  assert(argc > 2);

  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  my_size = malloc(world_size*sizeof(*my_size));
  displ = malloc(world_size*sizeof(*displ));
  displ[0] = 0;

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // printf("world sizehaha: %d", world_size);
  if (my_rank == 0){
    // printf("world sizehaha: %d", world_size);
    load(argv[1], &n, &vals);
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  rem = n % world_size;
  
  for (int i =0; i<world_size; i++){
    if (i < rem) my_size[i] = 1+n/world_size;
    else my_size[i] = n/world_size;
  }

  for(int i = 1; i < world_size; i++){
    displ[i] = displ[i-1] + my_size[i-1];
  }

  my_data = malloc(my_size[my_rank]*sizeof(*my_data));
  my_datas = calloc(world_size, sizeof(*my_datas));

  s = MPI_Wtime();

  MPI_Scatterv(
    vals,
    my_size,
    displ,
    MPI_UNSIGNED,
    my_data,
    my_size[my_rank],
    MPI_UNSIGNED,
    0,
    MPI_COMM_WORLD);

  if (my_rank == 0){
    free(vals);
  }

  sort_local(my_size[my_rank], my_data);

  // STEP 3
  // if(my_rank == 0){
  //   vals = malloc(n*sizeof(*vals));
  // }

  // MPI_Gatherv(
  //   my_data,
  //   my_size[my_rank],
  //   MPI_UNSIGNED,
  //   vals,
  //   my_size,
  //   displ,
  //   MPI_UNSIGNED,
  //   0,
  //   MPI_COMM_WORLD);
  vals = malloc(n*sizeof(*vals));
  int *pos_array = malloc(my_size[my_rank]*sizeof(*pos_array));
  MPI_Allgatherv(
    my_data,
    my_size[my_rank],
    MPI_UNSIGNED,
    vals,
    my_size,
    displ,
    MPI_UNSIGNED,
    MPI_COMM_WORLD);


  for (int i = displ[my_rank]; i < displ[my_rank]+my_size[my_rank]; i++){
    int sum = 0;
    for (int j = 0; j < world_size; j++){
      // if(my_rank == 0){
      //   printf("aaa: %d\n", sum);
      // }
      sum = sum + ryan_binary_search(displ[j], displ[j]+my_size[j]-1, vals, vals[i]);
    }
    // if(my_rank == 0){
    //   printf("-------------");
    // }
    pos_array[i-displ[my_rank]] = sum;
  }
  
  int *pos_arrays;
  if(my_rank == 0){
    pos_arrays = malloc(n*sizeof(*pos_arrays));
  }

  MPI_Gatherv(
    pos_array,
    my_size[my_rank],
    MPI_INT,
    pos_arrays,
    my_size,
    displ,
    MPI_INT,
    0,
    MPI_COMM_WORLD);
  
  if(my_rank == 0){
    unsigned int *result = calloc(n, sizeof(*result));
    for (int i = 0; i<n;i++){
      result[pos_arrays[i]] = vals[i];
      // if(result[i]==0)
      //   result[pos_arrays[i]] = vals[i];
      // else{
      //   int prev = i+1;
      //   while(1){
      //     if(result[prev] == 0) break;
      //     else prev = prev+1;
      //   }
      //   result[prev] = vals[i];
      // }
    }
    e = MPI_Wtime();
    print_time(e-s);
    print_numbers(argv[2], n, result);
  }

  // end step3


  // free(my_data);
  MPI_Finalize();
  return EXIT_SUCCESS;
}
