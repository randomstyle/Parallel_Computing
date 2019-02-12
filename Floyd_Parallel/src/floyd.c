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

#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
/* define access directions for matrices */
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]
#define a_res(R,C) a_res[ROWMJR(R,C,ln,n)]
 
#define INF 99999.0
 
// Solves the all-pairs shortest path problem using Floyd Warshall algorithm
void floydWarshall (int const n, float * a, float ** const res_p) {
    /* dist[][] will be the output matrix that will finally have the shortest 
      distances between every pair of vertices */
    int i, j, k;

    // float * a_res = malloc(n*n*sizeof(*a_res));
    // assert(a_res);
    /* Initialize the solution matrix same as input graph matrix. Or 
       we can say the initial values of shortest distances are based
       on shortest paths considering no intermediate vertex. */
    // for (i = 0; i < n; i++) {
    //     for (j = 0; j < n; j++) {
    //         a_res(i,j) = a(i,j);
    //     }
    // }
    /* Add all vertices one by one to the set of intermediate vertices.
      ---> Before start of a iteration, we have shortest distances between all
      pairs of vertices such that the shortest distances consider only the
      vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
      ----> After the end of a iteration, vertex no. k is added to the set of
      intermediate vertices and the set becomes {0, 1, 2, .. k} */
    for (k = 0; k < n; k++) {
        // Pick all vertices as source one by one
        for (i = 0; i < n; i++) {
            // Pick all vertices as destination for the
            // above picked source
            for (j = 0; j < n; j++) {
                // If vertex k is on the shortest path from
                // i to j, then update the value of dist[i][j]
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
    clock_t ts, te;

    if(argc < 3){
        printf("Invalid number of arguments.\nUsage: dijkstra <graph> <output_file>.\n");
        return EXIT_FAILURE;
    }

    load(argv[1], &n, &a);

    ts = clock();
    floydWarshall(n, a, &res);
    te = clock();

    print_numbers(argv[2],n,res);
    printf("Operation Time: %0.04fs\n", (double)(te-ts)/CLOCKS_PER_SEC);

    return 0;
}



// // C Program for Floyd Warshall Algorithm
// #include<stdio.h>
 
// // Number of vertices in the graph
// #define V 4
 
// /* Define Infinite as a large enough value. This value will be used
//   for vertices not connected to each other */
// #define INF 99999
 
// // A function to print the solution matrix
// void printSolution(int dist[][V]);
 
// // Solves the all-pairs shortest path problem using Floyd Warshall algorithm
// void floydWarshall (int graph[][V])
// {
//     /* dist[][] will be the output matrix that will finally have the shortest 
//       distances between every pair of vertices */
//     int dist[V][V], i, j, k;
 
//     /* Initialize the solution matrix same as input graph matrix. Or 
//        we can say the initial values of shortest distances are based
//        on shortest paths considering no intermediate vertex. */
//     for (i = 0; i < V; i++)
//         for (j = 0; j < V; j++)
//             dist[i][j] = graph[i][j];
 
//     /* Add all vertices one by one to the set of intermediate vertices.
//       ---> Before start of a iteration, we have shortest distances between all
//       pairs of vertices such that the shortest distances consider only the
//       vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
//       ----> After the end of a iteration, vertex no. k is added to the set of
//       intermediate vertices and the set becomes {0, 1, 2, .. k} */
//     for (k = 0; k < V; k++)
//     {
//         // Pick all vertices as source one by one
//         for (i = 0; i < V; i++)
//         {
//             // Pick all vertices as destination for the
//             // above picked source
//             for (j = 0; j < V; j++)
//             {
//                 // If vertex k is on the shortest path from
//                 // i to j, then update the value of dist[i][j]
//                 if (dist[i][k] + dist[k][j] < dist[i][j])
//                     dist[i][j] = dist[i][k] + dist[k][j];
//             }
//         }
//     }
 
//     // Print the shortest distance matrix
//     printSolution(dist);
// }
 
// /* A utility function to print solution */
// void printSolution(int dist[][V])
// {
//     printf ("Following matrix shows the shortest distances"
//             " between every pair of vertices \n");
//     for (int i = 0; i < V; i++)
//     {
//         for (int j = 0; j < V; j++)
//         {
//             if (dist[i][j] == INF)
//                 printf("%7s", "INF");
//             else
//                 printf ("%7d", dist[i][j]);
//         }
//         printf("\n");
//     }
// }
 
// // driver program to test above function
// int main()
// {
//     /* Let us create the following weighted graph
//             10
//        (0)------->(3)
//         |         /|\
//       5 |          |
//         |          | 1
//        \|/         |
//        (1)------->(2)
//             3           */
//     int graph[V][V] = { {0, 5, 3, 10},
//                         {5, 0, 3, 4},
//                         {3, 3, 0, 1},
//                         {10, 4, 1, 0}
//                       };
 
//     // Print the solution
//     floydWarshall(graph);
//     return 0;
// }