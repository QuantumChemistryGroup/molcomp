#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include "mc_types.h"
#include "mc_io.h"

static int mc_isBlank(char *line) {
    char *ch;
    int is_blank = -1;
    for (ch = line; *ch != '\0'; ++ch) {
        if (!isspace((unsigned char)*ch)) {
            is_blank = 0;
            break;
        }
    }
    return is_blank;
}

static char **MC_LINES(FILE *fp, long int *N) {
    char **list = malloc(1 * sizeof(char *));
    list[0] = malloc(100 * sizeof(char));
    char str[300];
    char *ptr;
    long int IT = 0;
    bool Skip2line = false;

    while (fgets(str, sizeof(str), fp) != NULL) {
        bool Skip = false;
        if ((Skip2line == true) && (IT == 1)) {
            Skip = true;
        }
        if (Skip == false) {
            if (mc_isBlank(str) == 0) {
                ptr = strtok(str, " \n\t\r");
                while (ptr != NULL) {
                    snprintf(list[*N], 100, "%s", ptr);
                    if (IT == 0 && (*N) == 0) {
                        if (atoi(ptr) != 0) {
                            Skip2line = true;
                        }
                    }
                    (*N)++;
                    list = realloc(list, (size_t)((*N) + 1) * sizeof(char *));
                    list[*N] = malloc(100 * sizeof(char));
                    ptr = strtok(NULL, " \n\t\r");
                }
            } else {
                strcpy(list[*N], "\n");
                (*N)++;
                list = realloc(list, (size_t)((*N) + 1) * sizeof(char *));
                list[*N] = malloc(100 * sizeof(char));
            }
        }
        IT = IT + 1;
    }
    return list;
}

void MC_FREAD(FILE *fp, McAtom **coordinates, long int *NATOMS, bool VER) {
    rewind(fp);
    long int N = 0;
    char **list = MC_LINES(fp, &N);

    long int i;
    long int NA = 0;
    long int j = 0;
    for (i = 0; i < N; i++) {
        char s[300];
        strcpy(s, list[i]);
        if (i == 0) {
            NA = atoi(s);
            if (NA == 0) {
                if (VER) printf("Warning: NATOMS in your xyz at first line is either 0 or missing.\n");
            } else {
                i = 1;
            }
        }
        if (i >= 0) {
            if (j == 0) {
                if (isalpha((unsigned char)list[i][0])) {
                    memset(&(*coordinates)[*NATOMS], 0, sizeof(McAtom));
snprintf((*coordinates)[*NATOMS].El, sizeof((*coordinates)[*NATOMS].El), "%s", list[i]);
                    j++;
                } else {
                    j = 0;
                }
            } else if (j > 0 && j < 4) {
                if (isalpha((unsigned char)list[i][0])) {
                    printf("Something is wrong with your coordinates.\n");
                    printf("Stop the code\n");
                    printf("%s\n", list[i]);
                    exit(0);
                } else {
                    (*coordinates)[*NATOMS].xyz[j - 1] = atof(list[i]);
                    j++;
                }
            }
            if (j == 4) {
                (*NATOMS)++;
                (*coordinates) = realloc((*coordinates), (size_t)((*NATOMS) + 1) * sizeof(McAtom));
                j = 0;
            }
        }
    }
    if (NA != (*NATOMS)) {
        if (NA == 0) {
            if (VER) printf("Warning: No number of atoms in the first line.\n");
        } else {
            printf("NA = %li NATOMS = %li\n", NA, (*NATOMS));
            printf("Stop the code\n");
            exit(0);
        }
    }
    for (i = 0; i < N + 1; i++) {
        free(list[i]);
    }
    free(list);
}

void MC_FPRINT(char *fileout, McAtom *coordinates, long int NATOMS, double RMSD) {
    FILE *fw = fopen(fileout, "w");
    if (!fw) return;
    fprintf(fw, " %li\n", NATOMS);
    fprintf(fw, " RMSD = %f\n", RMSD);
    for (long int i = 0; i < NATOMS; i++) {
        fprintf(fw, "%3s %12.8f %12.8f %12.8f\n",
                coordinates[i].El,
                coordinates[i].xyz[0],
                coordinates[i].xyz[1],
                coordinates[i].xyz[2]);
    }
    fclose(fw);
}

void MC_TXYZREAD(FILE *fp, McAtom **coordinates, long int NATOMS) {
    char str[300];
    char *ptr;
    long int i = 0, index = 0;

    if (fgets(str, 300, fp) != NULL) {}
    ptr = strtok(str, " ");
    ptr = strtok(NULL, " ");
    ptr = strtok(NULL, " ");
    char T[64] = "";
    if (ptr) {
        strncpy(T, ptr, sizeof(T) - 1);
    }

    if (fgets(str, 300, fp) != NULL) {}

    for (i = 0; i < NATOMS; i++) {
        (*coordinates)[i].nn = 0;
        ptr = strtok(str, " ");
        index = 0;
        while (ptr != NULL) {
            if (index == 5 && strcmp(T, "MM2") == 0) {
                (*coordinates)[i].t = atoi(ptr);
            } else {
                (*coordinates)[i].tu = NULL;
            }
            if (index > 5) {
            long int nb = strtol(ptr, NULL, 10);
            if (nb >= 1 && nb <= NATOMS && (*coordinates)[i].nn < 23) {
            (*coordinates)[i].n[(*coordinates)[i].nn] = nb;
            (*coordinates)[i].nn++;
            }
            }
            ptr = strtok(NULL, " ");
            index++;
        }
        if (fgets(str, 300, fp) != NULL) {}
    }
}
