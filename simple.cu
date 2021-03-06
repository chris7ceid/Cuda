#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>     
#include <sys/types.h>  
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>    
#include <errno.h>
#include <assert.h>

#define numThreads 1024

#define MAX_CHAR_PER_LINE 128

#define NONE 0
#define FIRST 1
#define LAST 2
#define BOTH 3 

double* file_read(char *filename,int  *numCoords,int  *numObjs) 
{
  double *objects;
  int     i, j, len;
    //ssize_t numBytesRead;
  int done=0; 
  FILE *infile;
  char *line, *ret;
  int   lineLen;

	//don't skip lines or attributes for this project
  int lines_to_skip=0; 
  int attr_to_skip=0;

  if ((infile = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Error: no such file (%s)\n", filename);
    return NULL;
  }

    /* first find the number of objects */
  lineLen = MAX_CHAR_PER_LINE;
  line = (char*) malloc(lineLen);
  assert(line != NULL);

  (*numCoords) = 0;

  while (fgets(line, lineLen, infile) != NULL) {
            /* check each line to find the max line length */
    while (strlen(line) == lineLen-1) {
                /* this line read is not complete */
      len = strlen(line);
      fseek(infile, -len, SEEK_CUR);

                /* increase lineLen */
      lineLen += MAX_CHAR_PER_LINE;
      line = (char*) realloc(line, lineLen);
      assert(line != NULL);

      ret = fgets(line, lineLen, infile);
      assert(ret != NULL);
    }

    if (strtok(line, " \t\n") != 0)
      (*numCoords)++;
  }

  (*numCoords)-=lines_to_skip;

  if((*numCoords)<=0)
  {
    fprintf(stderr, "Error: No objects found\n");
    return NULL;
  }

  rewind(infile);

	/*find the number of attributes*/  
  (*numObjs)=0;

  fgets(line, lineLen, infile);

  char * pch;
  pch=strtok(line, ",;");

  while (pch != NULL )
  {

    pch = strtok (NULL, ",;");
    (*numObjs)++;
  }

  if(attr_to_skip!=NONE)
  {
    (*numObjs)--;
    if(attr_to_skip==BOTH)
      (*numObjs)--;
  }

  rewind(infile);


    /* allocate space for objects and read all objects */
  len = (*numCoords) * (*numObjs);
  objects    = (double*)malloc( len * sizeof(double));
  assert(objects != NULL);



    /* read all objects */

  for(i=0;i<lines_to_skip;i++)
   fgets(line, lineLen, infile);

 i=0;
 j=0;

 while (fgets(line, lineLen, infile) != NULL) 
 {
   pch=strtok(line, ",;");
   while (pch != NULL && j<(*numObjs))
   {
    if(attr_to_skip%2==1 && j==0 && done==0)
    {
      done=1;
      pch = strtok (NULL, ",;");
      continue;                      
    }
    objects[i*(*numObjs)+j]=atof(pch);
    pch = strtok (NULL, ",;");
    j++;
  }
  i++;
  j=0;
  done=0;
}

assert(i == *numCoords);

fclose(infile);
free(line);


return objects;
}

__global__ void medval(double *M_d,double *w_d,int rows, int cols)
{
	int row = blockIdx.x*blockDim.x + threadIdx.x;
	double sum = 0.0;
	if (row < rows){
		for (int i=0;i<cols;i++){
			sum += M_d[row*cols + i];
		}
		w_d[row] = sum/cols;
	}	 
} 

__global__ void MatrixVecMul(double* M_d, double* x_d, double* z_d, double* w_d, int rows, int cols){
	int row = blockIdx.x*blockDim.x + threadIdx.x;
	double sum = 0.0;
	if (row < rows){	
		for (int i=0;i<cols;i++){
			sum += (M_d[row*cols + i] - w_d[row])*x_d[i];
		}
		z_d[row] = sum;
	}
} 



__global__ void TMatrixVecMul(double *M_d, double* x_d, double* z_d, double* k_d, double* w_d, int rows, int cols, float* norm_d){
	int row = blockIdx.x*blockDim.x + threadIdx.x;
	double sum = 0.0;
	if (row < cols){	
		for (int i=0;i<rows;i++){
			sum += (M_d[i*cols + row] - w_d[i])*z_d[i];
		}
		k_d[row] = sum;
		atomicAdd(norm_d,k_d[row]*k_d[row]);
	}

}

__global__ void DivideNorm(double* k_d, float* norm, int cols){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index < cols) k_d[index] = k_d[index]/sqrt(*norm);
}

__global__ void CalculateEps(double* k_d, double* x_d, double* e_d, float* eps, int cols){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index < cols){
   e_d[index] = k_d[index] - x_d[index];
   atomicAdd(eps,e_d[index]*e_d[index]);
 }
}


int main(int argc, char **argv){
	cudaSetDevice(0);
	int rows,cols;
  double* M_h = file_read(argv[1],&rows,&cols);
  int size = rows*cols*sizeof(double);
  double* x_h = (double*)malloc(cols*sizeof(double));
  double* z_h = (double*)malloc(rows*sizeof(double));
  double* k_h = (double*)malloc(cols*sizeof(double));
  float* e = (float*)malloc(sizeof(float));
  float* norm_h = (float*)malloc(sizeof(float));

  float* norm_d;
  float* eps;
  double* e_d;
  double* M_d;
  double* x_d;
  double* k_d;
  double* z_d;
  double* w_d;



  for (int i=0;i<cols;i++) x_h[i] = 1;


  cudaMalloc((void**) &M_d, size);
  cudaMalloc((void**) &x_d, cols*sizeof(double));
  cudaMalloc((void**) &z_d, rows*sizeof(double));
  cudaMalloc((void**) &k_d, cols*sizeof(double));
  cudaMalloc((void**) &e_d, cols*sizeof(double));
  cudaMalloc((void**) &w_d, rows*sizeof(double));
  cudaMalloc((void**) &norm_d, sizeof(float));
  cudaMalloc((void**) &eps, sizeof(float));

  cudaMemcpy(M_d, M_h, size, cudaMemcpyHostToDevice);
  cudaMemcpy(x_d, x_h, cols*sizeof(double), cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float elapsedTime;
  cudaEventRecord(start);
  medval<<<(rows - 1)/numThreads + 1,numThreads>>>(M_d,w_d,rows,cols);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  for (int i=0;i<100;i++){
    *norm_h = 0;
    *e = 0;
    cudaMemcpy(norm_d, norm_h, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(eps, e, sizeof(float), cudaMemcpyHostToDevice);
    if (i !=0){
     for (int k=0;k<cols;k++) x_h[k] = k_h[k];
       cudaMemcpy(x_d, x_h, cols*sizeof(double), cudaMemcpyHostToDevice);
   }
   cudaEventRecord(start);
   MatrixVecMul<<<(rows - 1)/numThreads + 1,numThreads>>>(M_d,x_d,z_d,w_d,rows,cols);

   TMatrixVecMul<<<(cols - 1)/numThreads + 1,numThreads>>>(M_d,x_d,z_d,k_d, w_d, rows,cols,norm_d);

   DivideNorm<<<(cols - 1)/numThreads + 1,numThreads>>>(k_d, norm_d, cols);

   CalculateEps<<<(cols - 1)/numThreads + 1,numThreads>>>(k_d, x_d, e_d, eps, cols);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   float milliseconds = 0;
   cudaEventElapsedTime(&milliseconds, start, stop);
   elapsedTime += milliseconds;

   cudaMemcpy(k_h,k_d,cols*sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(e, eps, sizeof(float), cudaMemcpyDeviceToHost);

   if (sqrt(*e) < 0.000001) break;
 }
 printf("Elapsed Time: %f ms\n", elapsedTime);
 FILE *f = fopen("results_simple.csv", "w");

 for (int i=0;i<cols;i++) fprintf(f,"%.7f%s", k_h[i],(i<cols-1)?",":"");
  fprintf(f,"\n");
fclose(f);

cudaFree(M_d);
cudaFree(x_d);
cudaFree(z_d);
cudaFree(k_d);
cudaFree(w_d);
cudaFree(norm_d);
cudaFree(e_d);
cudaFree(eps);

free(M_h);
free(x_h);
free(z_h);
free(k_h);
free(norm_h);
free(e);
return 0;
}

