#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>     
#include <sys/types.h>  
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>    
#include <errno.h>
#include <assert.h>

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

int main(int argc, char **argv){
	int rows,cols;
	double* M = file_read(argv[1],&rows,&cols);

	double* w = (double*)malloc(rows*sizeof(double));
	double* x = (double*)malloc(cols*sizeof(double));
	double* z = (double*)malloc(rows*sizeof(double));
	double* k = (double*)malloc(cols*sizeof(double));

	double elapsedTime = 0;
	int i,j,l;
	for (i=0;i<cols;i++) x[i] = 1;

	clock_t begin = clock();
	for (i = 0; i < rows; ++i)
	{
		w[i] = 0;
		for (l=0;l<cols;l++){
			w[i] += M[i*cols + l];
		}
		w[i] = w[i]/cols;	
	}
	clock_t end = clock();

	elapsedTime += (double)(end - begin);


	for (i=0;i<100;i++){
		begin = clock();
		for (j=0;j<rows;j++){
			z[j] = 0;
			for (l=0;l<cols;l++){
				z[j] += (M[j*cols + l] - w[j])*x[l];
			}
		}

		for (j=0;j<cols;j++){
			k[j] = 0;
			for (l=0;l<rows;l++){
				k[j] += (M[l*cols + j] - w[l])*z[l];
			}
		}

		double norm = 0;
		for (j=0;j<cols;j++){
			norm += k[j]*k[j];
		}
		norm = sqrt(norm);

		for (j=0;j<cols;j++){
			k[j] = k[j]/norm;
		}

		double e = 0;
		for (j=0;j<cols;j++){
			e += (k[j] - x[j])*(k[j] - x[j]);
		}

		if (sqrt(e) < 0.000001) break;

		end = clock();
		elapsedTime += (double)(end - begin);

		for (j=0;j<cols;j++){
			x[j] = k[j];
		}
	}

	printf("Elapsed Time: %f ms\n", elapsedTime/CLOCKS_PER_SEC*1000);
	FILE *f = fopen("results_seq.csv", "w");
	for (i=0;i<cols;i++) fprintf(f,"%.7f%s", k[i],(i<cols-1)?",":"");
	fprintf(f,"\n");
	fclose(f);


	return 0;
}