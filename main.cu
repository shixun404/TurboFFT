#include "include/turbofft.h"
// #include "include/code_gen/generated/fft_radix_2_logN_8_upload_0.cuh"
#include "include/code_gen/generated/fft_radix_2_logN_2_upload_0.cuh"
// #define DataType nv_bfloat162
#define DataType double2
// #define DataType float2
// #define DataType half2



void test_turbofft(DataType* input_d, DataType* output_d, DataType* output_turbofft, 
                    std::vector<long long int> param, long long int bs, int ntest){
    long long int N = (1 << param[0]), bs, Ni, WorkerFFTSize;
    long long int shared_size[3], griddims[3], blockdims[3]; 
    int kernel_launch_times = param[1];
    float gflops, elapsed_time;
    cudaEvent_t fft_begin, fft_end;
    
    for(int i = 0; i < kernel_launch_times; ++i){
        bs = param[5 + i];
        Ni = (1 << param[2 + i]); 
        WorkerFFTSize = param[8 + i]; 
        shared_size[i] = Ni * bs * sizeof(DataType);
        griddims[i] = N / (Ni * bs);
        blockdims[i] = (Ni * bs) / WorkerFFTSize;
        cudaFuncSetAttribute(turboFFTArr[logN][i], cudaFuncAttributeMaxDynamicSharedMemorySize, shared_size[i]);
    }
    
    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);

    cudaEventRecord(fft_begin);
    for (int j = 0; j < ntest; ++j){
        for(int i = 0; i < kernel_launch_times; ++i){
            turboFFTArr[logN][i]<<<griddims[i], blockdims[i], shared_size[i]>>>(input_d, output_d);
        }
        cudaDeviceSynchronize();
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);

    
    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    
    printf("turboFFT finished: T=%8.3fms, FLOPS=%8.3fGFLOPS\n", elapsed_time, gflops);
    
    checkCudaErrors(cudaMemcpy((void*)output_turbofft, (void*)output_d, N * sizeof(DataType), cudaMemcpyDeviceToHost));
}


int main(int argc, char *argv[]){
    if (argc != 3) {
        std::cerr << "Usage: program_name N bs" << std::endl;
        return 1;
    }

    long long N = std::atoi(argv[1]); // Convert first argument to integer
    long long bs = std::atoi(argv[2]); // Convert second argument to integer
    
    DataType* input, *output_turbofft, *output_cufft;
    DataType* input_d, *output_d;
    int ntest = 10;

    std::vector<std::vector<int>> params;
    params = utils::load_parameters();


    // Verification
    utils::initializeData<DataType>(input, input_d, output_d, output_turbofft, output_cufft, N, bs + 3);
    test_turbofft(input_d, output_d, output_turbofft, params[N], bs, 1);
    profiler::cufft::test_cufft<DataType>(input_d, output_d, output_cufft, N, bs, 1);
    
    utils::compareData<DataType>(output_turbofft, output_cufft, N, 1e-5);

    // Profiling
    test_turbofft(input_d, output_d, output_turbofft, params[N], bs, ntest);
    profiler::cufft::test_cufft<DataType>(input_d, output_d, output_cufft, N, bs, ntest);
    
    cudaFree(input_d);
    cudaFree(output_d);
    return 0;
}

// profiler::cufft::test_cufft_ft<DataType>(input_d, output_d, output_cufft, input_d + N * (bs + 2),
//                                          input_d + N * (bs + 1), output_d + N * (bs + 2),   N, bs, ntest);
// profiler::cufft::test_cufft_ft<DataType>(input_d, output_d, output_cufft, input_d + N * (bs + 2),
//                                          input_d + N * (bs + 1), output_d + N * (bs + 2),   N, bs, ntest);
