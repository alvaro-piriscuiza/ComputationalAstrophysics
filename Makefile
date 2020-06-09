COMPILER	= gcc
FLAGS		= -O0 -fopenmp
LIBRARIES	= -lm

project00: kepler.c
	   $(COMPILER) $(FLAGS) kepler.c -o kepler
					 
project01: leapfrog.c
           $(COMPILER) $(FLAGS) leapfrog.c -o leapfrog
					 
project02: mandelbrot.c
           $(COMPILER) $(FLAGS) mandelbrot.c -o mandelbrot
					 
project03: leapfrog.c
           $(COMPILER) $(FLAGS) rk4.c -o rk4
					 
					 
clean:
	rk4 rm project00 project01 project02 project03 
