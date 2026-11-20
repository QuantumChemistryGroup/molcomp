// (Minimum) Assignment Problem by Hungarian Algorithm
// taken from Knuth's Stanford Graphbase
// https://uta.edu

#include <stdio.h>
#include <stdlib.h>
#include "mc_hungarian.h"

#define INF (0x7FFFFFFF)
#define verbose (0)

void MC_HUNGMAIN(long long int size1, long long int size2, long long int COST[][size2], long int **MATCH) {
    long long int (*Array)[size2] = COST;
    char (*Result)[size2] = malloc(size1 * sizeof(*Result)); // used as boolean

    long int i;
    long int y, x;

    mc_hungarian(size1, size2, Array, Result);
    i = 0;
    for (y = 0; y < size1; ++y) {
        for (x = 0; x < size2; ++x) {
            if (Result[y][x]) {
                (*MATCH)[i] = x;
                i++;
            }
        }
    }
    free(Result);
}

void mc_hungarian(long long int size1, long long int size2, long long int Array[][size2], char Result[][size2]) {
    long long int i, j;
    int False = 0;
    int True = 1;

    long long int m = size1, n = size2;
    long long int k;
    long long int l;
    long long int s;

    long long int *col_mate = malloc(size1 * sizeof(long long int));
    long long int *unchosen_row = malloc(size1 * sizeof(long long int));
    long long int *row_dec = malloc(size1 * sizeof(long long int));
    for (i = 0; i < size1; i++) {
        col_mate[i] = 0;
        unchosen_row[i] = 0;
        row_dec[i] = 0;
    }

    long long int *row_mate = malloc(size2 * sizeof(long long int));
    long long int *parent_row = malloc(size2 * sizeof(long long int));
    long long int *col_inc = malloc(size2 * sizeof(long long int));
    long long int *slack = malloc(size2 * sizeof(long long int));
    long long int *slack_row = malloc(size2 * sizeof(long long int));

    for (i = 0; i < size2; i++) {
        row_mate[i] = 0;
        parent_row[i] = 0;
        col_inc[i] = 0;
        slack[i] = 0;
        slack_row[i] = 0;
    }

    long long int t;
    long long int q;
    long long int unmatched;
    long long int cost = 0;

    for (i = 0; i < size1; ++i)
        for (j = 0; j < size2; ++j)
            Result[i][j] = False;

    // Begin subtract column minima in order to start with lots of zeroes
    for (l = 0; l < n; l++) {
        s = Array[0][l];
        for (k = 1; k < n; k++)
            if (Array[k][l] < s)
                s = Array[k][l];
        cost += s;
        if (s != 0)
            for (k = 0; k < n; k++)
                Array[k][l] -= s;
    }

    // Begin initial state
    t = 0;
    for (l = 0; l < n; l++) {
        row_mate[l] = -1;
        parent_row[l] = -1;
        col_inc[l] = 0;
        slack[l] = INF;
    }
    for (k = 0; k < m; k++) {
        s = Array[k][0];
        for (l = 1; l < n; l++)
            if (Array[k][l] < s)
                s = Array[k][l];
        row_dec[k] = s;
        for (l = 0; l < n; l++)
            if (s == Array[k][l] && row_mate[l] < 0) {
                col_mate[k] = l;
                row_mate[l] = k;
                if (verbose)
                    printf("matching col %lld==row %lld\n", l, k);
                goto row_done;
            }
        col_mate[k] = -1;
        if (verbose)
            printf("node %lld: unmatched row %lld\n", t, k);
        unchosen_row[t++] = k;
row_done:
        ;
    }

    // Begin Hungarian algorithm
    if (t == 0)
        goto done;
    unmatched = t;
    while (1) {
        if (verbose)
            printf("Matched %lld rows.\n", m - t);
        q = 0;
        while (1) {
            while (q < t) {
                k = unchosen_row[q];
                s = row_dec[k];
                for (l = 0; l < n; l++)
                    if (slack[l]) {
                        long long int del;
                        del = Array[k][l] - s + col_inc[l];
                        if (del < slack[l]) {
                            if (del == 0) {
                                if (row_mate[l] < 0)
                                    goto breakthru;
                                slack[l] = 0;
                                parent_row[l] = k;
                                if (verbose)
                                    printf("node %lld: row %lld==col %lld--row %lld\n",
                                           t, row_mate[l], l, k);
                                unchosen_row[t++] = row_mate[l];
                            } else {
                                slack[l] = del;
                                slack_row[l] = k;
                            }
                        }
                    }
                q++;
            }

            // Begin introduce a new zero into the matrix
            s = INF;
            for (l = 0; l < n; l++)
                if (slack[l] && slack[l] < s)
                    s = slack[l];
            for (q = 0; q < t; q++)
                row_dec[unchosen_row[q]] += s;
            for (l = 0; l < n; l++)
                if (slack[l]) {
                    slack[l] -= s;
                    if (slack[l] == 0) {
                        k = slack_row[l];
                        if (verbose)
                            printf("Decreasing uncovered elements by %lld produces zero at [%lld,%lld]\n",
                                   s, k, l);
                        if (row_mate[l] < 0) {
                            for (j = l + 1; j < n; j++)
                                if (slack[j] == 0)
                                    col_inc[j] += s;
                            goto breakthru;
                        } else {
                            parent_row[l] = k;
                            if (verbose)
                                printf("node %lld: row %lld==col %lld--row %lld\n", t, row_mate[l], l, k);
                            unchosen_row[t++] = row_mate[l];
                        }
                    }
                } else
                    col_inc[l] += s;
        }
breakthru:
        // Begin update the matching
        if (verbose)
            printf("Breakthrough at node %lld of %lld!\n", q, t);
        while (1) {
            j = col_mate[k];
            col_mate[k] = l;
            row_mate[l] = k;
            if (verbose)
                printf("rematching col %lld==row %lld\n", l, k);
            if (j < 0)
                break;
            k = parent_row[j];
            l = j;
        }
        if (--unmatched == 0)
            goto done;
        // Begin get ready for another stage
        t = 0;
        for (l = 0; l < n; l++) {
            parent_row[l] = -1;
            slack[l] = INF;
        }
        for (k = 0; k < m; k++)
            if (col_mate[k] < 0) {
                if (verbose)
                    printf("node %lld: unmatched row %lld\n", t, k);
                unchosen_row[t++] = k;
            }
    }
done:
    // Begin doublecheck the solution
    for (k = 0; k < m; k++)
        for (l = 0; l < n; l++)
            if (Array[k][l] < row_dec[k] - col_inc[l])
                exit(0);
    for (k = 0; k < m; k++) {
        l = col_mate[k];
        if (l < 0 || Array[k][l] != row_dec[k] - col_inc[l])
            exit(0);
    }
    k = 0;
    for (l = 0; l < n; l++)
        if (col_inc[l])
            k++;
    if (k > m)
        exit(0);

    for (i = 0; i < m; ++i) {
        Result[i][col_mate[i]] = True;
    }
    for (k = 0; k < m; ++k) {
        for (l = 0; l < n; ++l) {
            Array[k][l] = Array[k][l] - row_dec[k] + col_inc[l];
        }
    }
    for (i = 0; i < m; i++)
        cost += row_dec[i];
    for (i = 0; i < n; i++)
        cost -= col_inc[i];

    free(col_mate);
    free(unchosen_row);
    free(row_dec);
    free(row_mate);
    free(parent_row);
    free(col_inc);
    free(slack);
    free(slack_row);
}
