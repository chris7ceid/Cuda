all: sequential.c simple.cu optimized.cu
	gcc -o seq sequential.c -lm
	nvcc -o simple simple.cu -arch sm_20
	nvcc -o opt optimized.cu -arch sm_20
