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

#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
/* define access directions for matrices */
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]
#define a_temp(R,C) a_temp[ROWMJR(R,C,ln,n)]
#define ttt(R,C) ttt[ROWMJR(R,C,ln,n)]

static void
load(
  char const * const filename,
  int * const np,
  float ** const ap,
  int num_node_n_proc,
  int my_rank
)
{
  int i, j, n, ret;
  FILE * fp=NULL;
  float * a;
  float * a_temp;
  // int * heads;
  
  /* open the file */
  fp = fopen(filename, "r");
  assert(fp);

  /* get the number of nodes in the graph */
  ret = fscanf(fp, "%d", &n);
  assert(1 == ret);

  int temp;
  int rem = n%num_node_n_proc;

  if (my_rank >= num_node_n_proc-rem && rem != 0) temp = n/num_node_n_proc + 1;
  else temp = n/num_node_n_proc;
  
  a = calloc(n*temp, sizeof(*a));
  assert(a);

  if (my_rank == 0){
    int temp1;
    for (int k=0; k<num_node_n_proc; ++k){
      if (k >= num_node_n_proc-rem && rem != 0) temp1 = n/num_node_n_proc + 1;
      else temp1 = n/num_node_n_proc;
      // printf("%d", temp1);
      a_temp = calloc(n*temp1, sizeof(*a_temp));
      assert(a_temp);

      for (i=0; i<temp1; ++i) {
        for (j=0; j<n; ++j) {
          if (k == 0){
            ret = fscanf(fp, "%f", &a(i,j));
            assert(1 == ret);
            // printf("a is %f\n", a(i,j));
          }
          else{
            ret = fscanf(fp, "%f", &a_temp(i,j));
            assert(1 == ret);
            // printf("a_temp is %f\n", a_temp(i,j));
          }
        }
        
      }
      
      if (k != 0) MPI_Send(a_temp, n*temp1, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
      free(a_temp);
    }
  }
  else MPI_Recv(a, n*temp, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  /* close file */
  ret = fclose(fp);
  assert(!ret);

  /* record output values */
  *np = n;
  *ap = a;
}

static void
dijkstra(
  int const s,
  int const n,
  float const * const a,
  float ** const lp,
  int my_rank,
  int num_node_n_proc
)
{
  int i, j;
  struct float_int {
    float l;
    int u;
  } min;
  char * m;
  float * l;
  float gMin[2];
  float * all_min;
  float * ttt;
  float minArr[2];
  float * l_part;
  

  m = calloc(n, sizeof(*m));
  assert(m);

  l = malloc(n*sizeof(*l));
  assert(l);
  
  int tee = n/num_node_n_proc;
  // int rem = n%num_node_n_proc;
  // if (my_rank >= num_node_n_proc-rem && rem != 0) temp = n/num_node_n_proc + 1;
  // else temp = n/num_node_n_proc;
  l_part = calloc(tee, sizeof(*l_part));
  assert(l_part);


  all_min = calloc(2*num_node_n_proc,sizeof(*all_min));
  assert(all_min);

  if(my_rank == 0){
    for (i=0; i<n; ++i) {
      l[i] = a(s,i);
    }
  }

  MPI_Bcast(l, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  m[s] = 1;
  min.u = -1; /* avoid compiler warning */

  for (i=1; i<n; ++i) {
    min.l = INFINITY;
    // if(my_rank ==2 && i==90){
      
    //   for(int aa = 0;aa<n; aa++){
    //     printf("m is %d\n", m[aa]);
    //   }
    // }
    // if(my_rank == 3){
    //   for(int aa = 0;aa<n; aa++){
    //     printf("%d %s\n",i, m[aa]);
    //   }
    // }
    /* find local minimum */
    for (j=(n/num_node_n_proc)*my_rank; j<(n/num_node_n_proc)*(my_rank+1); ++j) {
      if (!m[j] && l[j] < min.l) {
        min.l = l[j];
        min.u = j;
      }
    }

    minArr[0] = min.l;
    minArr[1] = (float)min.u;
    
    MPI_Gather(minArr, 2, MPI_FLOAT, all_min, 2, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    // if (my_rank==0){
    //   for(int cc =0; cc<num_node_n_proc*2;cc+=2){
    //     printf("%d %f %f\n",i, minArr[cc], minArr[cc+1]);
    //   }
    // }

    // if(my_rank == 0) {
    //   printf("%d\n", i);
    //   for(int ii=0; ii<20;ii+=2){
    //     printf("%f, %f: \n",all_min[ii],all_min[ii+1]);
    //   }
    // }

    // if(my_rank ==2 && i==1){
    //   for(int aa = 0;aa<n; aa++){
    //     printf("%d %f %f\n",i, all_min[aa], all_min[aa+1]);
    //   }
    // }
    if(my_rank == 0){
      gMin[0] = INFINITY;
      gMin[1] = INFINITY;
      for(int ii=0; ii<2*num_node_n_proc; ii+=2){
        if (all_min[ii] < gMin[0]){
          gMin[0] = all_min[ii];
          gMin[1] = all_min[ii+1];
          
        }
      }
    }

    MPI_Bcast(gMin, 2, MPI_FLOAT, 0, MPI_COMM_WORLD);

    min.l = gMin[0];
    min.u = (int)gMin[1];

    // if(my_rank==0) printf("loweest %d %d\n",i, min.u);

    m[min.u] = 1;
    // printf("m %f\n", min.l);



    // -----------//
    // if(i==1) {
    //   ttt = calloc(n*n, sizeof(*ttt));
    //   MPI_Gather(a, n*n/num_node_n_proc, MPI_FLOAT, ttt, n*n/num_node_n_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // }

    // if(my_rank==0){
    //   for (j=0; j<n; ++j) {
    //     if (!m[j] && min.l+ttt(j,min.u) < l[j]){
    //       l[j] = min.l+ttt(j,min.u);
    //     }
    //   }
    // }
    // MPI_Bcast(l, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // -----------//

    // for (j=0; j<n; ++j) {
    //   if (!m[j+(n/num_node_n_proc)*my_rank] && min.l + a(j+(n/num_node_n_proc)*my_rank, min.u) < l[j+(n/num_node_n_proc)*my_rank]){
    //     l_part[j] = min.l+a(j+(n/num_node_n_proc)*my_rank, min.u);
    //     // printf("%f\n", l_part[j-10*my_rank]);
    //     // if(ccc==0) printf("%d %f %f\n", my_rank, a(0,min.u), l[j] );
    //   }
    // }

    // MPI_Allgather(l_part, n/num_node_n_proc, MPI_FLOAT, l, n/num_node_n_proc, MPI_FLOAT, MPI_COMM_WORLD);
    // free(l_part);
    //copy l to l_part 000003000009
    

    int tempp = n/num_node_n_proc;
    //check all l ;
    //printf("%f\n", l[10]);
    for (j=tempp*my_rank; j<tempp*(my_rank+1); ++j) {
      
      l_part[j-tempp*my_rank] = l[j];
      
      int abc = j - (tempp*my_rank);
      // printf("%d %f\n", j, a(abc, min.u));
      if (!m[j] && min.l + a(abc, min.u) < l[j]){
        // if(j==10){
        //   for(i=0;i<10;i++){
        //     printf("%f %f\n", a(i, min.u), min.l);
        //   }
        // }
        
        l_part[j-tempp*my_rank] = min.l + a(abc, min.u);
        // printf("%f\n", l_part[j-tempp*my_rank]);
        // if(ccc==0) printf("%d %f %f\n", my_rank, a(0,min.u), l[j] );
      }
    }
    // for (j=tempp*my_rank; j<tempp*(my_rank+1); ++j) {
    //   printf("l part %d %f\n",j, l_part[j-tempp*my_rank]);
    // }
    MPI_Allgather(l_part, tempp, MPI_FLOAT, l, tempp, MPI_FLOAT, MPI_COMM_WORLD);
    // if(my_rank==0){
    //   for(i=0;i<n;i++){
    //     printf("l %d is %f \n",i,l[i]);
    //   }
    // }
    //free(l_part);

  }

  free(m);
  *lp = l;
}

static void
print_time(double const seconds)
{
  printf("Operation Time: %0.04fs\n", seconds);
}

static void
print_numbers(
  char const * const filename,
  int const n,
  float const * const numbers)
{
  int i;
  FILE * fout;

  /* open file */
  if(NULL == (fout = fopen(filename, "w"))) {
    fprintf(stderr, "error opening '%s'\n", filename);
    abort();
  }

  /* write numbers to fout */
  for(i=0; i<n; ++i) {
    fprintf(fout, "%10.4f\n", numbers[i]);
  }

  fclose(fout);
}

int
main(int argc, char ** argv)
{
  
  int n;
  double ts, te;
  float * a, * l;
  // int * heads;

  MPI_Init(NULL, NULL);
  int world_size;
  int my_rank;

  int num_node_n_proc = atoi(argv[4]) * atoi(argv[5]);

  if(argc < 4){
     printf("Invalid number of arguments.\nUsage: dijkstra <graph> <source> <output_file> <num_node> <num_proc>.\n");
     return EXIT_FAILURE;
  }


  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // my_size = malloc(world_size*sizeof(*my_size));
  // displ = malloc(world_size*sizeof(*displ));
  // displ[0] = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  

  load(argv[1], &n, &a, num_node_n_proc, my_rank);
  // if (my_rank == 0){
  //   printf("zero%f\n", a(0,0));
  // }
  // else{
  //   printf("%d %f\n",my_rank, a(0,0));
  // }

  ts = MPI_Wtime();
  dijkstra(atoi(argv[2]), n, a, &l, my_rank, num_node_n_proc);
  if(my_rank == 0){
    te = MPI_Wtime();
    print_time(te-ts);
    print_numbers(argv[3], n, l);
  }
  
  free(a);
  free(l);
  MPI_Finalize();
  return EXIT_SUCCESS;
}
