#include "include/turbofft.h"
// #include "include/code_gen/generated/fft_radix_2_logN_8_upload_0.cuh"
#include "include/code_gen/generated/fft_radix_2_logN_2_upload_0.cuh"
// #define DataType nv_bfloat162
#define DataType double2
// #define DataType float2
// #define DataType half2



void test_turbofft(DataType* input_d, DataType* output_d, DataType* output_turbofft, 
                    long long int N){
    dim3 gridDim(1, 1, 1); 
    dim3 blockDim(1, 1, 1);
    // turbofft::fft::thread::fft<DataType, turbofft::Tensor<DataType, 1, 2>><<<gridDim, blockDim>>>(input_d, output_d);
    // turbofft::fft::thread::fft<DataType><<<gridDim, blockDim>>>(input_d, output_d);
    long long int shared_size = N * sizeof(DataType);

    // cudaFuncSetAttribute(turbofft::fft::thread::fft<DataType, 8>, cudaFuncAttributeMaxDynamicSharedMemorySize, shared_size);
    cudaFuncSetAttribute(fft_radix_2_logN_2_dim_0, cudaFuncAttributeMaxDynamicSharedMemorySize, shared_size);
    
    
    // turbofft::fft::thread::fft<DataType, 8><<<gridDim, blockDim, shared_size>>>(input_d, output_d);
    // turbofft::fft::thread::fft<DataType><<<gridDim, blockDim>>>(input_d, output_d);
    // fft_radix_2_logN_8_dim_0<<<gridDim, blockDim>>>(input_d, output_d);
    fft_radix_2_logN_2_dim_0<<<gridDim, blockDim>>>(input_d, output_d);
    cudaDeviceSynchronize();
    printf("%d\n",  N * sizeof(DataType));
    checkCudaErrors(cudaMemcpy((void*)output_turbofft, (void*)output_d, N * sizeof(DataType), cudaMemcpyDeviceToHost));
}


int main(int argc, char *argv[]){
    DataType* input, *output_turbofft, *output_cufft;
    DataType* input_d, *output_d;
    long long int N = 1 << 2, bs = 1;
    int ntest = 1;

    if (argc < 2) bs = 1;
    else bs = std::atoi(argv[1]);
    printf("N=%d, bs=%d\n", N, bs);
    utils::initializeData<DataType>(input, input_d, output_d, output_turbofft, output_cufft, N, bs + 3);
    test_turbofft(input_d, output_d, output_turbofft, N);

    
    // profiler::cufft::test_cufft_ft<DataType>(input_d, output_d, output_cufft, input_d + N * (bs + 2),
    //                                          input_d + N * (bs + 1), output_d + N * (bs + 2),   N, bs, ntest);

    profiler::cufft::test_cufft<DataType>(input_d, output_d, output_cufft, N, bs, ntest);

    // profiler::cufft::test_cufft_ft<DataType>(input_d, output_d, output_cufft, input_d + N * (bs + 2),
    //                                          input_d + N * (bs + 1), output_d + N * (bs + 2),   N, bs, ntest);

    
    utils::compareData<DataType>(output_turbofft, output_cufft, N, 1e-5);
    // printData(output_turbofft, N);
    // printData(output_cufft, N);
    cudaFree(input_d);
    cudaFree(output_d);
    return 0;
}