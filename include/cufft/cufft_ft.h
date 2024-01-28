// #include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cublas_v2.h>  
#include <cuda_fp16.h>
#include <cuda_bf16.h>
#include "reduction.cuh"

namespace profiler{
namespace cufft{
template<typename DataType>
void test_cufft_ft(DataType* input_d, DataType* output_d, DataType* output_cufft,
                DataType* e_d, DataType* x_i_c, DataType* x_o_c, 
                long long int N, size_t bs, int ntest, int M);
template<typename DataType>
void test_cufft_ft_gemv(DataType* input_d, DataType* output_d, DataType* output_cufft,
                DataType* e_d, DataType* x_i_c, DataType* x_o_c, 
                long long int N, size_t bs, int ntest);

template<>
void test_cufft_ft<float2>(float2* input_d, float2* output_d, float2* output_cufft,
                            float2* e_d, float2* x_i_c, float2* x_o_c,
                            long long int N, size_t bs, int ntest, int M) {
    cufftHandle plan, plan_ft;
    cublasHandle_t handle;         
    float gflops, elapsed_time, mem_bandwidth;
    cuComplex alpha = {1, 1}, beta = {0, 0};
    cudaEvent_t fft_begin, fft_end;
    cublasCreate(&handle);  
    dim3 gridDim1((N + 255) / 256, bs / M, 1);

    checkCudaErrors(cufftCreate(&plan));
    checkCudaErrors(cufftCreate(&plan_ft));

    checkCudaErrors(cufftPlan1d(&plan, N, CUFFT_C2C, bs + 1));
    // checkCudaErrors(cufftPlan1d(&plan_ft, N, CUFFT_C2C, 16));
    // cublasSgemv(handle, CUBLAS_OP_N, N * 2, bs / 16, (float*)&(alpha), 
    //                                 reinterpret_cast<float*>(input_d), N * 2, 
    //                                 reinterpret_cast<float*>(e_d), 1, (float*)&(beta), 
    //                                 reinterpret_cast<float*>(x_i_c), 1);

    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int i = 0; i < ntest; ++i){
        // cublasSgemv(handle, CUBLAS_OP_N, N, bs, (float*)&(alpha), 
        //                             reinterpret_cast<float*>(input_d), N, 
        //                             reinterpret_cast<float*>(e_d), 1, (float*)&(beta), 
        //                              reinterpret_cast<float*>(x_i_c), 1);
        cublasSgemv(handle, CUBLAS_OP_T, N, bs, (float*)&(alpha), 
                                    reinterpret_cast<float*>(input_d), N, 
                                    reinterpret_cast<float*>(e_d), 1, (float*)&(beta), 
                                     reinterpret_cast<float*>(x_i_c), 1);
        
        my_checksum<<<gridDim1, 256>>>(N, M, reinterpret_cast<float*>(input_d),
                                             reinterpret_cast<float*>(x_i_c));
                
        checkCudaErrors(cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(input_d), 
                     reinterpret_cast<cufftComplex*>(output_d), 
                     CUFFT_FORWARD));
        // checkCudaErrors(cufftExecC2C(plan_ft, reinterpret_cast<cufftComplex*>(input_d), 
        //              reinterpret_cast<cufftComplex*>(output_d), 
        //              CUFFT_FORWARD));
        
      
        cublasSgemv(handle, CUBLAS_OP_T, N, bs, (float*)&(alpha), 
                                    reinterpret_cast<float*>(output_d), N, 
                                    reinterpret_cast<float*>(e_d), 1, (float*)&(beta), 
                                    reinterpret_cast<float*>(x_o_c), 1);
      
      
        cudaDeviceSynchronize();
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    mem_bandwidth = (float)(N * bs * 8 * 2) / (elapsed_time) * 1000.f / 1000000000.f;
    printf("cuFFT_FT_online_mychecksum, %d, %d, %8.3f, %8.3f, %8.3f\n",  (int)log2f(N),  (int)log2f(bs), elapsed_time, gflops, mem_bandwidth);

    // printf("cuFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);

    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(float2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}

template<>
void test_cufft_ft<double2>(double2* input_d, double2* output_d, double2* output_cufft,
                            double2* e_d, double2* x_i_c, double2* x_o_c,
                            long long int N, size_t bs, int ntest, int M) {
    cufftHandle plan;
    float gflops, elapsed_time, mem_bandwidth;
    cudaEvent_t fft_begin, fft_end;
    cublasHandle_t handle;         
    cuDoubleComplex alpha = {1, 1}, beta = {0, 0};
    
    dim3 gridDim1((N + 255) / 256, bs / M, 1);
    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftPlan1d(&plan, N, CUFFT_Z2Z, bs));

    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int i = 0; i < ntest; ++i){
        // cublasDgemv(handle, CUBLAS_OP_N, N, bs, (double*)&(alpha), 
        //                             reinterpret_cast<double*>(input_d), N, 
        //                             reinterpret_cast<double*>(e_d), 1, (double*)&(beta), 
        //                              reinterpret_cast<double*>(x_i_c), 1);
        // cudaDeviceSynchronize();
        cublasDgemv(handle, CUBLAS_OP_T, N, bs, (double*)&(alpha), 
                                    reinterpret_cast<double*>(input_d), N, 
                                    reinterpret_cast<double*>(e_d), 1, (double*)&(beta), 
                                     reinterpret_cast<double*>(x_i_c), 1);
        cudaDeviceSynchronize();
        my_checksum<<<gridDim1, 256>>>(N, M, reinterpret_cast<double*>(input_d),
                                             reinterpret_cast<double*>(x_i_c));
                
        cudaDeviceSynchronize();
        checkCudaErrors(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(input_d), 
                        reinterpret_cast<cufftDoubleComplex*>(output_d), 
                        CUFFT_FORWARD));
      cudaDeviceSynchronize();
                cublasDgemv(handle, CUBLAS_OP_T, N, bs, (double*)&(alpha), 
                                    reinterpret_cast<double*>(output_d), N, 
                                    reinterpret_cast<double*>(e_d), 1, (double*)&(beta), 
                                    reinterpret_cast<double*>(x_o_c), 1);
      cudaDeviceSynchronize();
      
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    // printf("cuFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);
    mem_bandwidth = (float)(N * bs * 16 * 2) / (elapsed_time) * 1000.f / 1000000000.f;
    printf("cuFFT, %d, %d, %8.3f, %8.3f, %8.3f\n",  (int)log2f(N),  (int)log2f(bs), elapsed_time, gflops, mem_bandwidth);



    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * bs * sizeof(double2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}

}
}