#!/bin/bash
set -e

CC=${CC:-gcc}
CFLAGS="-std=c11 -O3 -march=native -fopenmp -Wall -Wextra"

$CC -c mc_atominfo.c    $CFLAGS
$CC -c mc_inertia.c     $CFLAGS
$CC -c mc_kabshangle.c  $CFLAGS
$CC -c mc_getconnect.c  $CFLAGS
$CC -c mc_stdfunc.c     $CFLAGS
$CC -c mc_hungarian.c   $CFLAGS
$CC -c mc_renumb_hung.c $CFLAGS
$CC -c mc_core.c        $CFLAGS
$CC -c mc_io.c          $CFLAGS
$CC -c main.c           $CFLAGS
$CC -c moldedupl.c      $CFLAGS

OBJS="mc_atominfo.o mc_inertia.o mc_kabshangle.o mc_getconnect.o mc_stdfunc.o mc_hungarian.o mc_renumb_hung.o mc_core.o mc_io.o"

$CC -o molcomp   $OBJS main.o      -lm -O3 -march=native -fopenmp
$CC -o moldedupl $OBJS moldedupl.o -lm -O3 -march=native -fopenmp

rm -f *.o

echo "Built molcomp and moldedupl successfully."

