#include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cublas_v2.h>  
#include <cuda_fp16.h>
#include <cuda_bf16.h>


namespace profiler{
namespace cufft{
template<typename DataType>
void test_cufft_ft(DataType* input_d, DataType* output_d, DataType* output_cufft,
                DataType* e_d, DataType* x_i_c, DataType* x_o_c, 
                long long int N, size_t bs, int ntest);

template<>
void test_cufft_ft<float2>(float2* input_d, float2* output_d, float2* output_cufft,
                            float2* e_d, float2* x_i_c, float2* x_o_c,
                            long long int N, size_t bs, int ntest) {
    cufftHandle plan, plan_ft;
    cublasHandle_t handle;         
    float gflops, elapsed_time;
    cuComplex alpha = {1, 1}, beta = {0, 0};
    cudaEvent_t fft_begin, fft_end;
    cublasCreate(&handle);  
    

    checkCudaErrors(cufftCreate(&plan));
    checkCudaErrors(cufftCreate(&plan_ft));

    checkCudaErrors(cufftPlan1d(&plan, N, CUFFT_C2C, bs + 1));
    checkCudaErrors(cufftPlan1d(&plan_ft, N, CUFFT_C2C, 16));


    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int i = 0; i < ntest; ++i){
        // if(cublasCgemm(handle, CUBLAS_OP_N,CUBLAS_OP_N, N, 1, bs, &alpha, 
        //                             reinterpret_cast<cuComplex*>(input_d), N, 
        //                             reinterpret_cast<cuComplex*>(e_d), bs, &beta, 
        //                             reinterpret_cast<cuComplex*>(x_i_c), N) != CUBLAS_STATUS_SUCCESS){
        //     printf("cuBLASCGEMM ERROR!\n");
        //     return;
        // }
        // cublasCgemv(handle, CUBLAS_OP_N, N, bs / 8, &alpha, 
        //                             reinterpret_cast<cuComplex*>(input_d), N, 
        //                             reinterpret_cast<cuComplex*>(e_d), 1, &beta, 
        //                             reinterpret_cast<cuComplex*>(x_i_c), 1);
        // cublasSgemv(handle, CUBLAS_OP_N, N * 2, bs / 16, (float*)&(alpha), 
        //                             reinterpret_cast<float*>(input_d), N * 2, 
        //                             reinterpret_cast<float*>(e_d), 1, (float*)&(beta), 
        //                             reinterpret_cast<float*>(x_i_c), 1);
        checkCudaErrors(cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(input_d), 
                     reinterpret_cast<cufftComplex*>(output_d), 
                     CUFFT_FORWARD));
        checkCudaErrors(cufftExecC2C(plan_ft, reinterpret_cast<cufftComplex*>(input_d), 
                     reinterpret_cast<cufftComplex*>(output_d), 
                     CUFFT_FORWARD));
        
        // if(cublasCgemm(handle, CUBLAS_OP_N,CUBLAS_OP_N, N, 1, bs, &alpha, 
        //                             reinterpret_cast<cuComplex*>(output_d), N, 
        //                             reinterpret_cast<cuComplex*>(e_d), bs, &beta, 
        //                             reinterpret_cast<cuComplex*>(x_o_c), N) != CUBLAS_STATUS_SUCCESS){
        //     printf("cuBLASCGEMM ERROR!\n");
        //     return;
        // }
        // cublasCgemv(handle, CUBLAS_OP_N, N, bs / 8, &alpha, 
        //                             reinterpret_cast<cuComplex*>(output_d), N, 
        //                             reinterpret_cast<cuComplex*>(e_d), 1, &beta, 
        //                             reinterpret_cast<cuComplex*>(x_o_c), 1);
        // cublasSgemv(handle, CUBLAS_OP_N, N * 2, bs / 16, (float*)&(alpha), 
        //                             reinterpret_cast<float*>(output_d), N * 2, 
        //                             reinterpret_cast<float*>(e_d), 1, (float*)&(beta), 
        //                             reinterpret_cast<float*>(x_o_c), 1);
        // checkCudaErrors(cublasCgemm(handle, CUBLAS_OP_N,CUBLAS_OP_N, N, 1, bs, &alpha, output_d, N, e_d, bs, &beta, x_o_c, N));
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    printf("cuFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);

    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(float2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}


template<>
void test_cufft_ft<double2>(double2* input_d, double2* output_d, double2* output_cufft,
                            double2* e_d, double2* x_i_c, double2* x_o_c,
                            long long int N, size_t bs, int ntest) {
    cufftHandle plan;
    float gflops, elapsed_time;
    cudaEvent_t fft_begin, fft_end;
    

    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftPlan1d(&plan, N, CUFFT_Z2Z, bs));

    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int i = 0; i < ntest; ++i){
        checkCudaErrors(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(input_d), 
                        reinterpret_cast<cufftDoubleComplex*>(output_d), 
                        CUFFT_FORWARD));
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    printf("cuFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);


    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(double2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}

template<>
void test_cufft_ft<nv_bfloat162>(nv_bfloat162* input_d, nv_bfloat162* output_d, nv_bfloat162* output_cufft,
                                nv_bfloat162* e_d, nv_bfloat162* x_i_c, nv_bfloat162* x_o_c,
                                long long int N, size_t bs, int ntest) {
    cufftHandle plan;
    float gflops, elapsed_time;
    cudaEvent_t fft_begin, fft_end;
    size_t ws = 0;

    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftXtMakePlanMany(plan, 1, &N, NULL, 0, 0, CUDA_C_16BF, 
                    NULL, 0, 0, CUDA_C_16BF, bs, &ws, CUDA_C_16BF));


    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int i = 0; i < ntest; ++i){
        checkCudaErrors(cufftXtExec(plan, reinterpret_cast<nv_bfloat162*>(input_d), 
                        reinterpret_cast<nv_bfloat162*>(output_d), 
                        CUFFT_FORWARD));
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    printf("cuFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);



    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(nv_bfloat162), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}

template<>
void test_cufft_ft<half2>(half2* input_d, half2* output_d, half2* output_cufft,
                         half2* e_d, half2* x_i_c, half2* x_o_c,
                         long long int N, size_t bs, int ntest) {
    cufftHandle plan;
    float gflops, elapsed_time;
    cudaEvent_t fft_begin, fft_end;
    size_t ws = 0;

    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftXtMakePlanMany(plan, 1, &N, NULL, 0, 0, CUDA_C_16F, 
                    NULL, 0, 0, CUDA_C_16F, bs, &ws, CUDA_C_16F));

    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int i = 0; i < ntest; ++i){
        checkCudaErrors(cufftXtExec(plan, reinterpret_cast<half2*>(input_d), 
                        reinterpret_cast<half2*>(output_d), 
                        CUFFT_FORWARD));
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    printf("cuFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);


    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(half2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}
}
}