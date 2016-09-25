all:
	gcc -o cgSolver cgSolver.c -c99 -lm
clean:
	rm -rf *.o
