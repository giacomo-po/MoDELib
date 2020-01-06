
// /usr/local/cuda/bin/nvcc testSpeedUp.cu -o test -O3

// WARNING for OPTIMIZATION
// warning: compiling with nvcc -O3 filename.cu will pass the -O3 option to host code only.
// nvcc -Xptxas -O3,-v filename.cu
// https://stackoverflow.com/questions/43706755/how-can-i-get-the-nvcc-cuda-compiler-to-optimize-more

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <time.h>

#define N (1024*1024)
#define M (10000)
#define THREADS_PER_BLOCK 1024

void serial_add(double *a, double *b, double *c, int n, int m)
{
for(int index=0;index<n;index++)
{
for(int j=0;j<m;j++)
{
c[index] = a[index]*a[index] + b[index]*b[index];
}
}
}

__global__ void vector_add(double *a, double *b, double *c)
{
int index = blockIdx.x * blockDim.x + threadIdx.x;
for(int j=0;j<M;j++)
{
c[index] = a[index]*a[index] + b[index]*b[index];
}
}

int main()
{
clock_t start,end;

double *a, *b, *c;
int size = N * sizeof( double );

a = (double *)malloc( size );
b = (double *)malloc( size );
c = (double *)malloc( size );

for( int i = 0; i < N; i++ )
{
a[i] = b[i] = i;
c[i] = 0;
}

start = clock();
serial_add(a, b, c, N, M);

printf( "c[0] = %d\n",0,c[0] );
printf( "c[%d] = %d\n",N-1, c[N-1] );

end = clock();

float time1 = ((float)(end-start))/CLOCKS_PER_SEC;
printf("Serial: %f seconds\n",time1);

start = clock();
double *d_a, *d_b, *d_c;


cudaMalloc( (void **) &d_a, size );
cudaMalloc( (void **) &d_b, size );
cudaMalloc( (void **) &d_c, size );


cudaMemcpy( d_a, a, size, cudaMemcpyHostToDevice );
cudaMemcpy( d_b, b, size, cudaMemcpyHostToDevice );

vector_add<<< (N + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( d_a, d_b, d_c );

cudaMemcpy( c, d_c, size, cudaMemcpyDeviceToHost );


printf( "c[0] = %d\n",0,c[0] );
printf( "c[%d] = %d\n",N-1, c[N-1] );


free(a);
free(b);
free(c);
cudaFree( d_a );
cudaFree( d_b );
cudaFree( d_c );

end = clock();
float time2 = ((float)(end-start))/CLOCKS_PER_SEC;
printf("CUDA: %f seconds, Speedup: %f\n",time2, time1/time2);

return 0;
}
