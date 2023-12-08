#include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cuda_fp16.h>
#include <cuda_bf16.h>


namespace profiler{
namespace cufft{
template<typename DataType>
void test_cufft(DataType* input_d, DataType* output_d, 
                DataType* output_cufft, long long int N, size_t bs);

template<>
void test_cufft<float2>(float2* input_d, float2* output_d, 
                        float2* output_cufft, long long int N, size_t bs) {
    cufftHandle plan;

    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftPlan1d(&plan, N, CUFFT_C2C, bs));

    checkCudaErrors(cufftExecC2C(plan, reinterpret_cast<cufftComplex*>(input_d), 
                     reinterpret_cast<cufftComplex*>(output_d), 
                     CUFFT_FORWARD));

    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(float2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}
template<>
void test_cufft<double2>(double2* input_d, double2* output_d, 
                        double2* output_cufft, long long int N, size_t bs) {
    cufftHandle plan;

    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftPlan1d(&plan, N, CUFFT_Z2Z, bs));

    checkCudaErrors(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(input_d), 
                     reinterpret_cast<cufftDoubleComplex*>(output_d), 
                     CUFFT_FORWARD));

    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(double2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}

template<>
void test_cufft<nv_bfloat162>(nv_bfloat162* input_d, nv_bfloat162* output_d, 
                        nv_bfloat162* output_cufft, long long int N, size_t bs) {
    cufftHandle plan;
    size_t ws = 0;

    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftXtMakePlanMany(plan, 1, &N, NULL, 0, 0, CUDA_C_16BF, 
                    NULL, 0, 0, CUDA_C_16BF, bs, &ws, CUDA_C_16BF));

    checkCudaErrors(cufftXtExec(plan, reinterpret_cast<nv_bfloat162*>(input_d), 
                     reinterpret_cast<nv_bfloat162*>(output_d), 
                     CUFFT_FORWARD));

    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(nv_bfloat162), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}

template<>
void test_cufft<half2>(half2* input_d, half2* output_d, 
                        half2* output_cufft, long long int N, size_t bs) {
    cufftHandle plan;
    size_t ws = 0;

    checkCudaErrors(cufftCreate(&plan));

    checkCudaErrors(cufftXtMakePlanMany(plan, 1, &N, NULL, 0, 0, CUDA_C_16F, 
                    NULL, 0, 0, CUDA_C_16F, bs, &ws, CUDA_C_16F));

    checkCudaErrors(cufftXtExec(plan, reinterpret_cast<half2*>(input_d), 
                     reinterpret_cast<half2*>(output_d), 
                     CUFFT_FORWARD));

    checkCudaErrors(cudaMemcpy(output_cufft, output_d, N * sizeof(half2), 
                   cudaMemcpyDeviceToHost));

    checkCudaErrors(cufftDestroy(plan));
}
}
}