// #include "include/turbofft/tensor.h"
#include "include/turbofft/utils.h"
#include "include/turbofft/fft/thread/fft.h"
#include <stdio.h>
#include <cuda_runtime.h> 
#include <cufftXt.h>
#define DataType float2

void test_cufft(DataType* input_d, DataType* output_d, DataType* output_cufft, size_t N){
    cufftHandle plan;  
    cufftCreate(&plan);
    int res = cufftPlan1d(&plan, N, CUFFT_C2C, 1);
    res = cufftExecC2C(plan, (cufftComplex *)input_d, (cufftComplex *)output_d, CUFFT_FORWARD);
    cudaMemcpy((void*)output_cufft, (void*)output_d, N * sizeof(DataType), cudaMemcpyDeviceToHost);
    cufftDestroy(plan);
}

void test_turbofft(DataType* input_d, DataType* output_d, DataType* output_turbofft, 
                    size_t N){
    dim3 gridDim(1, 1, 1); 
    dim3 blockDim(1, 1, 1);
    // turbofft::fft::thread::fft<DataType, turbofft::Tensor<DataType, 1, 2>><<<gridDim, blockDim>>>(input_d, output_d);
    // turbofft::fft::thread::fft<DataType><<<gridDim, blockDim>>>(input_d, output_d);
    turbofft::fft::thread::fft<<<gridDim, blockDim>>>(input_d, output_d);
    cudaDeviceSynchronize();
    printf("%d\n",  N * sizeof(DataType));
    CUDA_CALLER(cudaMemcpy((void*)output_turbofft, (void*)output_d, N * sizeof(DataType), cudaMemcpyDeviceToHost));

}

void compareData(DataType* res, DataType *res_ref, size_t N, double error_bound, 
                bool printInfo=false){
    double rel_error = 0.;
    for(int i = 0; i < N; ++i){
        rel_error = abs((res[i].x - res_ref[i].x) / res_ref[i].x);
        if(rel_error > error_bound){
            printf("Error detected: res[%d].x=%f, res_ref[%d].x=%f, rel_error=%f\n", 
            i, res[i].x, i, res_ref[i].x, rel_error);
        }
        rel_error = abs((res[i].y - res_ref[i].y) / res_ref[i].y);
        if(rel_error > error_bound){
            printf("Error detected: res[%d].y=%f, res_ref[%d].y=%f, rel_error=%f\n", 
            i, res[i].y, i, res_ref[i].y, rel_error);
        }
    }
}
void printData(DataType* res, size_t N){
    double rel_error = 0.;
    for(int i = 0; i < N; ++i){
        printf("res[%d] = %f + %f j\n", i, res[i].x, res[i].y);
    }
}

void initializeData(DataType *&input, DataType *&input_d, DataType *&output_d, 
                    DataType *&output_turbofft, DataType *&output_cufft, size_t N){
    input = (DataType*)calloc(N, sizeof(DataType));
    output_turbofft = (DataType*)calloc(N, sizeof(DataType));
    output_cufft = (DataType*)calloc(N, sizeof(DataType));
    CUDA_CALLER(cudaMalloc((void**)&input_d, sizeof(DataType) * N));
    CUDA_CALLER(cudaMalloc((void**)&output_d, sizeof(DataType) * N));
    
    for(int i = 0; i < N; ++i){
        input[i].x = 1;
        input[i].y = 1;
    }

    CUDA_CALLER(cudaMemcpy((void*)input_d, (void*)input, N * sizeof(DataType), cudaMemcpyHostToDevice));
}

int main(){
    DataType* input, *output_turbofft, *output_cufft;
    DataType* input_d, *output_d;
    int N = 2;

    initializeData(input, input_d, output_d, output_turbofft, output_cufft, N);
    test_turbofft(input_d, output_d, output_turbofft, N);

    test_cufft(input_d, output_d, output_cufft, N);

    compareData(output_turbofft, output_cufft, N, 1e-5, false);
    printData(output_turbofft, N);
    printData(output_cufft, N);

    return 0;
}