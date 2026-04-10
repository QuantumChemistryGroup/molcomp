#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// this code is from
// https://www.sanfoundry.com/c-program-implement-fisher-yates-algorithm-array-shuffling/ 

 

static long int rand_int(long int n) {
    long int limit = RAND_MAX - RAND_MAX % n;
    long int rnd;
    do {
        rnd = rand();
    } 
    while (rnd >= limit);
    return rnd % n;
}

void shuffle(long int *array, long int n) {
    int i, j, tmp;
    for (i = n - 1; i > 0; i--) {
        j = rand_int(i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
   }
}
