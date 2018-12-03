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

__global__ void median(double *M_d,double *w_d,int rows, int cols)
{
	__shared__ double w_ds[numThreads];
	int row = blockIdx.x;

	int loop = (cols - 1)/numThreads + 1;
	double sum = 0.0;
	for (int k=0;k<loop;k++){ //Μερικά Αθροίσματα
		if ( row < rows && (numThreads*k + threadIdx.x) < cols ){
			sum += M_d[row*cols + numThreads*k + threadIdx.x];
		}
	}
	w_ds[threadIdx.x] = sum;
	__syncthreads();

	for (unsigned int k=blockDim.x/2; k>0; k>>=1) { //Reduction
		if (threadIdx.x < k) {
			w_ds[threadIdx.x] += w_ds[threadIdx.x + k];
		}
		__syncthreads();
	}
	w_d[row] = w_ds[0]/cols;

} 

__global__ void MatrixVecFirst(double* M_d, double* x_d, double* z_d, double* w_d, int rows, int cols){
	__shared__ double Pvalue[numThreads];
	__shared__ double w_ds;

	int row = blockIdx.x;
	if (!threadIdx.x) w_ds = w_d[row];
	__syncthreads();

	int loop = (cols - 1)/numThreads + 1;
	double sum = 0.0;
	for (int k=0;k<loop;k++){ //Μερικά Αθροίσματα
		if ( row < rows && (numThreads*k + threadIdx.x) < cols ) 
		{
			sum+= (M_d[row*cols + numThreads*k + threadIdx.x] - w_ds)*x_d[numThreads*k + threadIdx.x];

		}         
	}
	Pvalue[threadIdx.x]=sum;
	__syncthreads();

for (unsigned int k=blockDim.x/2; k>0; k>>=1) { //Reduction
		if (threadIdx.x < k) {
			Pvalue[threadIdx.x] += Pvalue[threadIdx.x + k];
		}
		__syncthreads();
	}
	z_d[row] = Pvalue[0]; //Save result
} 



__global__ void MatrixVecSecond(double *M_d, double* x_d, double* z_d, double* k_d, double* w_d, int rows, int cols){
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int Row = bx*numThreads + tx;
	double sum = 0.0;
	if (Row < cols){ //Κάθε thread μια γραμμή
		for (int k=0;k<rows;k++){ //rows αντί για cols λόγω ανάστροφου
			sum += (M_d[k*cols + Row] - w_d[k])*z_d[k]; //διαφορετικός μο
		}
		k_d[Row] = sum;
	}	
}

__global__ void CalculateNorm(double* k_d, float* norm, int cols){
	__shared__ double k_ds[numThreads];

	int tx = threadIdx.x;
	int index = blockIdx.x*numThreads + tx;

	if (index < cols) k_ds[tx] = k_d[index]*k_d[index];
	else k_ds[tx] = 0;
	__syncthreads();

	for (unsigned int k=blockDim.x/2; k>0; k>>=1) {
		if (tx < k) {
			k_ds[tx] += k_ds[tx + k];
		}
		__syncthreads();
	}

	if (tx == 0) atomicAdd(norm,k_ds[tx]);
}

__global__ void DivideByNorm(double* k_d, float* norm, int cols){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index < cols) k_d[index] = k_d[index]/sqrt(*norm);
}

__global__ void CalculateEps(double* k_d, double* x_d, double* e_d, int cols){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index < cols) e_d[index] = k_d[index] - x_d[index];
}


int main(int argc, char **argv){
	int rows,cols;

	double* M_h = file_read(argv[1],&rows,&cols);
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


	//Arxikopoiisi tou x_0
	for (int i=0;i<cols;i++) x_h[i] = 1;

	cudaSetDevice(0);
	cudaMalloc((void**) &M_d, rows*cols*sizeof(double));
	cudaMalloc((void**) &x_d, cols*sizeof(double));
	cudaMalloc((void**) &z_d, rows*sizeof(double));
	cudaMalloc((void**) &k_d, cols*sizeof(double));
	cudaMalloc((void**) &e_d, cols*sizeof(double));
	cudaMalloc((void**) &w_d, rows*sizeof(double));
	cudaMalloc((void**) &norm_d, sizeof(float));
	cudaMalloc((void**) &eps, sizeof(float));

	cudaMemcpy(M_d, M_h, rows*cols*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(x_d, x_h, cols*sizeof(double), cudaMemcpyHostToDevice);

	float elapsedTime;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start);
	median<<<rows,numThreads>>>(M_d,w_d,rows,cols);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&elapsedTime, start, stop);
	for (int i=0;i<100;i++){ //Endeiktiki timi to 100 gia periptoseis mh sigklisis
		*norm_h = 0;
		*e = 0;
		cudaMemcpy(norm_d, norm_h, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(eps, e, sizeof(float), cudaMemcpyHostToDevice);

		cudaEventRecord(start);
		MatrixVecFirst<<<rows,numThreads>>>(M_d,x_d,z_d,w_d,rows,cols); //z_d = (M_d - w_d*e')*x_d

		MatrixVecSecond<<<(cols + numThreads -1)/numThreads,numThreads>>>(M_d,x_d,z_d,k_d, w_d, rows,cols); //k_d = (M_d - w_d*e')'*z_d

		CalculateNorm<<<(cols + numThreads-1)/numThreads,numThreads>>>(k_d, norm_d, cols); //norm(k_d)

		DivideByNorm<<<(cols + numThreads -1)/numThreads,numThreads>>>(k_d, norm_d, cols); //x_k+1 = k_d/norm(k_d)

		CalculateEps<<<(cols + numThreads -1)/numThreads,numThreads>>>(k_d, x_d, e_d, cols); //x_k+1 - x_k

		CalculateNorm<<<(cols + numThreads -1)/numThreads,numThreads>>>(e_d, eps, cols); // Ypologismos eps = ||x_k+1 - x_k||
		cudaEventRecord(stop);
		cudaEventSynchronize(stop);

		float milliseconds = 0;
		cudaEventElapsedTime(&milliseconds, start, stop);
		elapsedTime += milliseconds;


		//Elegxos sigklisis
		cudaMemcpy(e, eps, sizeof(float), cudaMemcpyDeviceToHost);
		if (sqrt(*e) < 0.000001){
			cudaMemcpy(k_h,k_d,cols*sizeof(double),cudaMemcpyDeviceToHost);
			break;
		}

		//An den sigklinei, arxikopoiisi tou algorithmou me x_k+1
		cudaMemcpy(x_d,k_d,cols*sizeof(double),cudaMemcpyDeviceToDevice);
	}
	

	printf("Elapsed time: %f ms\n", elapsedTime);

	//Print to file
	FILE *f = fopen("results_opt.csv", "w");
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
