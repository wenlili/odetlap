COMPILER = nvcc --compiler-bindir /usr/bin/g++-4.9 -O2 -std=c++11 -D_FORCE_INLINES -arch=sm_35

all: main.o approximator.o solver.o
	$(COMPILER) -o program main.o approximator.o solver.o

main.o: main.cu approximator.h
	$(COMPILER) -c main.cu

approximator.o: approximator.cu approximator.h solver.h
	$(COMPILER) -c approximator.cu

solver.o: solver.cu solver.h
	$(COMPILER) -c solver.cu

run: program
	./program input/0001.txt input/nzm0001.txt 160 256 256 output/m0001.txt output/v0001.txt output/a0001.txt > output/o0001.txt &

test: count.py test.py
	python count.py output/m0001.txt
	python test.py input/0001.txt output/a0001.txt

runch: program
	./program input/201501ch4.txt input/nzm201501ch4.txt 24 180 360 output/m201501ch4.txt output/v201501ch4.txt output/a201501ch4.txt > output/o201501ch4.txt &

testch: count.py testnz.py
	python count.py output/m201501ch4.txt
	python testnz.py input/201501ch4.txt input/nzm201501ch4.txt output/a201501ch4.txt

runco: program
	./program input/201501co.txt input/nzm201501co.txt 24 180 360 output/m201501co.txt output/v201501co.txt output/a201501co.txt > output/o201501co.txt &

testco: count.py testnz.py
	python count.py output/m201501co.txt
	python testnz.py input/201501co.txt input/nzm201501co.txt output/a201501co.txt

rung: program
	./program input/201501gph.txt input/nzm201501gph.txt 24 180 360 output/m201501gph.txt output/v201501gph.txt output/a201501gph.txt > output/o201501gph.txt &

testg: count.py testnz.py
	python count.py output/m201501gph.txt
	python testnz.py input/201501gph.txt input/nzm201501gph.txt output/a201501gph.txt

runo: program
	./program input/201501o3.txt input/nzm201501o3.txt 24 180 360 output/m201501o3.txt output/v201501o3.txt output/a201501o3.txt > output/o201501o3.txt &

testo: count.py testnz.py
	python count.py output/m201501o3.txt
	python testnz.py input/201501o3.txt input/nzm201501o3.txt output/a201501o3.txt

runt: program
	./program input/201501temp.txt input/nzm201501temp.txt 24 180 360 output/m201501temp.txt output/v201501temp.txt output/a201501temp.txt > output/o201501temp.txt &

testt: count.py testnz.py
	python count.py output/m201501temp.txt
	python testnz.py input/201501temp.txt input/nzm201501temp.txt output/a201501temp.txt

sync:
	rsync makefile *.h *.cu *.py liw9@geoxeon.ecse.rpi.edu:3D1
