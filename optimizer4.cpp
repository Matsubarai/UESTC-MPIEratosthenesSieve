#include "mpi.h"
#include <math.h>
#include <stdio.h>

#define BLOCK_LOW(id, p, n) ((long long)(id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define ODD_TO_INDEX(value) (((value)-3)>>1)
#define INDEX_TO_ODD(index) (((index)<<1)+3)

//实验没说不给用pragma里的Ofast和加速指令！
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")

int main(int argc, char *argv[]) {
    int count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    int first;        /* Index of first multiple */
    int global_count; /* Global prime count */
    int i, j;
    int id;           /* Process ID number */
    int index;        /* Index of current prime */
    int low_value;    /* Lowest value on this proc */
    char *marked;     /* Portion of 2,...,'n' */
    int n;            /* Sieving from 2, ..., 'n' */
    int p;            /* Number of processes */
    int prime;        /* Current prime */
    int size;         /* Elements in 'marked' */
    int m;            /* Size of search list */
    int offset;
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

    n = atoi(argv[1]);
    m = (n - 1) / 2;

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */

    low_value = INDEX_TO_ODD(BLOCK_LOW(id, p, m));
    size = BLOCK_SIZE(id, p, m);

    /* Allocate this process's share of the array. */

    marked = (char *) calloc(size, sizeof(char));
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    primes_size = ODD_TO_INDEX((int)sqrt(n)) + 1;
    primes = (char *) calloc(primes_size, sizeof(char));
    if (primes == NULL) {
        printf("Cannot allocate enough memory\n");
        free(marked);
        MPI_Finalize();
        exit(1);
    }

    /* preprocess primes in 3..sqrt(n) */
    index = 0;
    prime = 3;
    do {
        for (i = ODD_TO_INDEX(prime * prime); i < primes_size; i += prime)
            primes[i] = 1;
        while (primes[++index]);
        prime = INDEX_TO_ODD(index);
    } while (prime * prime <= sqrt(n));

    for(i = 0; i < size; i += chunk){
        index = 0;
        prime = 3;
        low_value_chunk = INDEX_TO_ODD(ODD_TO_INDEX(low_value) + i);
            do {
                if (prime * prime > low_value_chunk)
                    first = ODD_TO_INDEX(prime * prime) - ODD_TO_INDEX(low_value_chunk);
                else {
                    offset = low_value_chunk % prime;   // temp value
                    if (!offset) first = 0;
                    else {
                        first = prime - offset;
                        if (!((low_value_chunk + first) & 1))
                            first += prime;
                        first >>= 1;
                    }
                }
                for (j = first + i; j < first + i + chunk && j < size; j += prime) marked[j] = 1;
                while (primes[++index]);
                prime = INDEX_TO_ODD(index);
            } while (prime * prime <= n);
    }

    count = size;
    for (i = 0; i < size; i++)
        count -= marked[i]; // delete 'if'
    if(p > 1) MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    else count = global_count;

    /* Stop the timer */

    elapsed_time += MPI_Wtime();

    /* Print the results */

    if (!id) {
        printf("There are %d primes less than or equal to %d\n",
               global_count + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}
