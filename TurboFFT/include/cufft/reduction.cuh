__global__ void my_checksum(int N, int M, float *a,  float *c)
{	
	a += (threadIdx.x + blockIdx.x * blockDim.x) + blockIdx.y * M * N;
	float result = 0;
	float tmp = *a;
	#pragma unroll
	for(int i = 0; i < M; ++i){
		result += tmp;
		a += N;
		tmp = *a;
	}
	c[threadIdx.x + blockIdx.x * blockDim.x + blockIdx.y * N] = result;
}

__global__ void my_checksum(int N, int M, double *a,  double *c)
{	
	a += (threadIdx.x + blockIdx.x * blockDim.x) + blockIdx.y * M * N;
	double result = 0;
	double tmp = *a;
	#pragma unroll
	for(int i = 0; i < M; ++i){
		result += tmp;
		a += N;
		tmp = *a;
	}
	c[threadIdx.x + blockIdx.x * blockDim.x + blockIdx.y * N] = result;
}