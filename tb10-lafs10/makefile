all:
	@ gcc -o cgSolver_v1 cgSolver_v1.c -std=c99 -lm -g -O3 -mavx -march=native -DLIKWID_PERFORM -L/home/soft/likwid/lib -I/home/soft/likwid/include -llikwid
	@ gcc -o cgSolver_v2 cgSolver_v2.c -std=c99 -lm -g -O3 -mavx -march=native -DLIKWID_PERFORM -L/home/soft/likwid/lib -I/home/soft/likwid/include -llikwid
clean:
	@ rm -rf *.o
