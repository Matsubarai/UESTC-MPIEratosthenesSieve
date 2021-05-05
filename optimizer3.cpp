#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b) ((a)<(b) ? (a) : (b))
#define BLOCK_LOW(id, p, n) ((long long)(id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index, p, n) (((p)*(index)+1)-1)/(n))
#define VALUE_TO_INDEX(value) (((value)-2)
#define INDEX_TO_VALUE(index) ((index)+2)
#define ODD_TO_INDEX(value) (((value)-3)/2)
#define INDEX_TO_ODD(index) (2*(index)+3)

int main(int argc, char *argv[]) {
    int count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    int first;        /* Index of first multiple */
    int global_count; /* Global prime count */
    int high_value;   /* Highest value on this proc */
    int i, j;
    int id;           /* Process ID number */
    int index;        /* Index of current prime */
    int low_value;    /* Lowest value on this proc */
    char *marked;     /* Portion of 2,...,'n' */
    int n;            /* Sieving from 2, ..., 'n' */
    int p;            /* Number of processes */
    int proc0_size;   /* Size of proc 0's subarray */
    int prime;        /* Current prime */
    int size;         /* Elements in 'marked' */
    int m;            /* Size of search list */
    char *primes;     /* Preprocessed primes */
    int primes_size;  /* Elements in 'primes' */
    int chunk;        /* chunk size for *marked to adapt cache size */
    int low_value_chunk; /* Lowest value in a chunk */

    freopen("/dev/null","w",stderr);
    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    chunk = (1 << 20) / p;  // 256KBX4 L2 Cache=1MB
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoi(argv[1]);
    m = ODD_TO_INDEX(n) + 1;

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */

    low_value = INDEX_TO_ODD(BLOCK_LOW(id, p, m));
    high_value = INDEX_TO_ODD(BLOCK_HIGH(id, p, m));
    size = BLOCK_SIZE(id, p, m);

    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = m / p;

    if (INDEX_TO_ODD(proc0_size - 1) < (int) sqrt((double) n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */

    marked = (char *) malloc(size);
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    for (i = 0; i < size; i++) marked[i] = 0;

    primes_size = ODD_TO_INDEX(sqrt(n)) + 1;
    primes = (char *) malloc(primes_size);
    if (primes == NULL) {
        printf("Cannot allocate enough memory\n");
        free(marked);
        MPI_Finalize();
        exit(1);
    }
    for (i = 0; i < primes_size; i++) primes[i] = 0;

    /* preprocess primes in 3..sqrt(n) */
    index = 0;
    prime = 3;
    do {
        for (i = ODD_TO_INDEX(prime * prime); i < primes_size; i += prime)
            primes[i] = 1;
        while (primes[++index]);
        prime = INDEX_TO_ODD(index);
    } while (prime * prime <= sqrt(n));

    for(i = 0; i < size; i += chunk){   // chunking
        index = 0;
        prime = 3;
        low_value_chunk = INDEX_TO_ODD(ODD_TO_INDEX(low_value) + i);
            do {
                if (prime * prime > low_value_chunk)
                    first = ODD_TO_INDEX(prime * prime) - ODD_TO_INDEX(low_value_chunk);
                else {
                    if (!(low_value_chunk % prime)) first = 0;
                    else {
                        first = prime - (low_value_chunk % prime);
                        if (!((low_value_chunk + first) & 1))
                            first += prime;
                        first >>= 1;
                    }
                }
                for (j = first + i; j < first + i + chunk && j < size; j += prime)  // update in first+i..min(first+i+chunk-1, size-1)
                    marked[j] = 1;
                while (primes[++index]);
                prime = INDEX_TO_ODD(index);
            } while (prime * prime <= n);
    }

    count = 0;
    for (i = 0; i < size; i++) {
        if (!marked[i])
            count++;
    }
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();

    /* Print the results */

    if (!id) {
        printf("There are %d primes less than or equal to %d\n", global_count + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}
