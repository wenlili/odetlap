COMPILER = nvcc -O2 -std=c++11 -D_FORCE_INLINES -arch=sm_35

all: main.o approximator.o solver.o
	$(COMPILER) -o program main.o approximator.o solver.o

main.o: main.cu approximator.h
	$(COMPILER) -c main.cu

approximator.o: approximator.cu approximator.h solver.h
	$(COMPILER) -c approximator.cu

solver.o: solver.cu solver.h
	$(COMPILER) -c solver.cu
