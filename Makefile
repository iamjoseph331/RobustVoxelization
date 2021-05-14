ALL: clean1 main
clean1:
	rm -rf *out*
	rm -rf *.gch 
	rm -rf prof
main:  
	g++ -std=c++11 -I . -O3 voxelize.cpp -o exe.out 
omp:
	g++ -std=c++11 -I . -O3 voxelize.cpp -o exe.out -fopenmp	
cuda:
	nvcc voxelize.cu -I . -O3 -o exe.out -Xcompiler -fopenmp -w 
debug: 
	g++ -std=c++11 -I . -g voxelize.cpp -o exe.out -pg
	gdb -q ./exe.out
ex:
		./exe.out < in > out
clean:
	rm -rf *out
prof:
	nvprof --unified-memory-profiling off ./exe.out -i out1.msh -o out -r 256 -p
