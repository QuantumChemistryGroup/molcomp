gcc -c atominfo.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only 
gcc -c fread.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only
gcc -c inertia.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -lgsl -lgslcblas
gcc -c kabsh.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -lgsl -lgslcblas
gcc -c kabshdist.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -lgsl -lgslcblas
gcc -c molcomp.c -fopenmp -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -lgsl -lgslcblas 
gcc -c renumb.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only
gcc -c hungarian.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only
gcc -c renumb-hung.c -fopenmp -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only
gcc -c shuffle.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only
gcc -c stdfunc.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only
gcc -c getconnect.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only
gcc -c txyzread.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only

gcc -o a.out molcomp.o fread.o atominfo.o inertia.o kabsh.o kabshdist.o renumb.o hungarian.o renumb-hung.o shuffle.o stdfunc.o getconnect.o txyzread.o -fopenmp -lgsl -lgslcblas -std=c11 -pedantic -lm -Ofast -fno-finite-math-only -Wall -static -static-libgcc # -fsanitize=address -g -fno-omit-frame-pointer -static-libasan #-static #-g

