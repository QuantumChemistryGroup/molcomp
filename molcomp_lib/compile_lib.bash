gcc -c atominfo.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c fread.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c inertia.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c kabsh.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c kabshdist.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c kabshangle.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c molcomp.c -fopenmp -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c renumb.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c hungarian.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c renumb-hung.c -fopenmp -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c shuffle.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c stdfunc.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c getconnect.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC
gcc -c txyzread.c -std=c11 -Wall -pedantic -Ofast -fno-finite-math-only -fPIC

gcc -shared -o libmolcomp.so molcomp.o fread.o atominfo.o inertia.o kabsh.o kabshdist.o kabshangle.o renumb.o hungarian.o renumb-hung.o shuffle.o stdfunc.o getconnect.o txyzread.o -fopenmp -lgsl -lgslcblas -lm -Ofast

