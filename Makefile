ALL: clean1 main
clean1:
	rm -rf *out
	rm -rf *.gch 
main:  
	g++ -std=c++11 -I . -O3 voxelize.cpp -o exe.out 
omp:
	g++ -std=c++11 -I . -O3 voxelize.cpp -o exe.out -fopenmp	
debug: 
	g++ -std=c++11 -I . -g voxelize.cpp -o exe.out -pg
	gdb -q ./exe.out
ex:
		./exe.out < in > out
clean:
	rm -rf *out