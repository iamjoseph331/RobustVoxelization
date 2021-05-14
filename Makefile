ALL: clean1 main
clean1:
	rm -rf *out
	rm -rf *.gch 
main:  
	g++ -std=c++11 -I . -O3 voxelize.cpp -o exe.out 
omp:
	g++ -std=c++11 -I . -O3 voxelize.cpp -o exe.out -fopenmp	
cuda:
	nvcc voxelize.cu -I . -o exe.out -Xcompiler -fopenmp -w 
debug: 
	g++ -std=c++11 -I . -g voxelize.cpp -o exe.out -pg
	gdb -q ./exe.out
ex:
		./exe.out < in > out
clean:
	rm -rf *out
cudatest:
	nvcc add.cu -I . -o cudatest.out -w
prof:
	nvprof --unified-memory-profiling off ./cudatest.out -i out1.msh -o out -r 128 
