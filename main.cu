#include "include/turbofft.h"
// #define DataType nv_bfloat162
// #define DataType double2
#define DataType double2
// #define DataType half2



void test_turbofft(DataType* input_d, DataType* output_d, DataType* output_turbofft, 
                    long long int N){
    dim3 gridDim(1, 1, 1); 
    dim3 blockDim(1, 32, 1);
    // turbofft::fft::thread::fft<DataType, turbofft::Tensor<DataType, 1, 2>><<<gridDim, blockDim>>>(input_d, output_d);
    // turbofft::fft::thread::fft<DataType><<<gridDim, blockDim>>>(input_d, output_d);
    long long int shared_size = (N * sizeof(DataType) / 16) * 17;

    cudaFuncSetAttribute(turbofft::fft::thread::fft<DataType, 8>, cudaFuncAttributeMaxDynamicSharedMemorySize, shared_size);
    turbofft::fft::thread::fft<DataType, 8><<<gridDim, blockDim, shared_size>>>(input_d, output_d);
    // turbofft::fft::thread::fft<DataType><<<gridDim, blockDim>>>(input_d, output_d);
    cudaDeviceSynchronize();
    printf("%d\n",  N * sizeof(DataType));
    checkCudaErrors(cudaMemcpy((void*)output_turbofft, (void*)output_d, N * sizeof(DataType), cudaMemcpyDeviceToHost));

}


int main(){
    DataType* input, *output_turbofft, *output_cufft;
    DataType* input_d, *output_d;
    int N = 256, bs=1;

    utils::initializeData<DataType>(input, input_d, output_d, output_turbofft, output_cufft, N);
    test_turbofft(input_d, output_d, output_turbofft, N);

    profiler::cufft::test_cufft<DataType>(input_d, output_d, output_cufft, N, bs);

    utils::compareData<DataType>(output_turbofft, output_cufft, N, 1e-5);
    // printData(output_turbofft, N);
    // printData(output_cufft, N);

    return 0;
}