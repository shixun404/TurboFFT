#include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cuda_fp16.h>
#include <cuda_bf16.h>


namespace profiler{
namespace cufft{
template<typename DataType>
void test_cufft(DataType* input_d, DataType* output_d, 
                DataType* output_cufft, long long int N, size_t bs, int ntest);

template<>
void test_cufft<float2>(float2* input_d, float2* output_d, 
                        float2* output_cufft, long long int N, size_t bs, int ntest) {
    cufftHandle plan;
    float gflops, elapsed_time, mem_bandwidth;
    cudaEvent_t fft_begin, fft_end;


    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftPlan1d(&plan, N, CUFFT_C2C, bs));

    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int i = 0; i < ntest; ++i){
        checkCudaErrors(cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(input_d), 
                        reinterpret_cast<cufftComplex*>(output_d), 
                        CUFFT_FORWARD));
        cudaDeviceSynchronize();
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    // printf("cuFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);
    mem_bandwidth = (float)(N * bs * 8 * 2) / (elapsed_time) * 1000.f / 1000000000.f;
    printf("cuFFT, %d, %d, %8.3f, %8.3f, %8.3f\n",  (int)log2f(N),  (int)log2f(bs), elapsed_time, gflops, mem_bandwidth);

    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * bs * sizeof(float2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}


template<>
void test_cufft<double2>(double2* input_d, double2* output_d, 
                        double2* output_cufft, long long int N, size_t bs, int ntest) {
    cufftHandle plan;
    float gflops, elapsed_time, mem_bandwidth;
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

template<>
void test_cufft<nv_bfloat162>(nv_bfloat162* input_d, nv_bfloat162* output_d, 
                        nv_bfloat162* output_cufft, long long int N, size_t bs, int ntest) {
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
void test_cufft<half2>(half2* input_d, half2* output_d, 
                        half2* output_cufft, long long int N, size_t bs, int ntest) {
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