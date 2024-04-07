
#include "include/TurboFFT.h"
    
template <typename DataType, int if_ft, int if_err>
void test_turbofft( DataType* input_d, DataType* output_d, DataType* output_turbofft,
                    DataType* twiddle_d, DataType* checksum, std::vector<long long int> param, 
                    long long int bs, int thread_bs, int ntest){
    long long int N = (1 << param[0]), threadblock_bs, Ni, WorkerFFTSize;
    long long int logN = param[0];
    long long int shared_size[3], griddims[3], blockdims[3]; 
    DataType* inputs[3] = {input_d, output_d, output_d + N * bs};
    DataType* outputs[3] = {output_d, output_d + N * bs, output_d};
    int kernel_launch_times = param[1];
    float gflops, elapsed_time, mem_bandwidth;
    cudaEvent_t fft_begin, fft_end;
    
    cublasHandle_t handle;      
    TurboFFT_Kernel_Entry<DataType, if_ft, if_err> entry;
    int M = 16;
    dim3 gridDim1((N + 255) / 256, bs / M, 1);
    for(int i = 0; i < kernel_launch_times; ++i){
        threadblock_bs = param[5 + i];
        Ni = (1 << param[2 + i]); 
        WorkerFFTSize = param[8 + i]; 
        shared_size[i] = Ni * threadblock_bs * sizeof(DataType);
        
        blockdims[i] = (Ni * threadblock_bs) / WorkerFFTSize;
        long long int shared_per_SM = 160 * 1024;
        shared_per_SM = 128 * 1024;
        griddims[i] = min(108 * min((2048 / blockdims[i]), (shared_per_SM / shared_size[i])), 
                ((N * bs) + (Ni * threadblock_bs) - 1) / (Ni * threadblock_bs));
        
        griddims[i] = ((((N * bs) + (Ni * threadblock_bs) - 1) / (Ni * threadblock_bs))) / thread_bs;
    
        cudaFuncAttributes attr;
        if(cudaFuncSetAttribute(entry.turboFFTArr[logN][i], cudaFuncAttributeMaxDynamicSharedMemorySize, shared_size[i]))
        printf("Set DynamicSharedMem failed\n");
        if(cudaFuncSetAttribute(entry.turboFFTArr[logN][i], cudaFuncAttributePreferredSharedMemoryCarveout, (shared_per_SM * 100) / (164 * 1024)))
        printf("Set smemCarveout failed\n");
        cudaError_t get_attr_res = cudaFuncGetAttributes (&attr, entry.turboFFTArr[logN][i] );
        if(get_attr_res != 0)
        printf("get_attr_res = %d\n", get_attr_res);
    }
    
    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);
    #pragma unroll
    for(int i = 0; i < kernel_launch_times; ++i){
        entry.turboFFTArr[logN][i]<<<griddims[i], blockdims[i], shared_size[i]>>>(inputs[i], outputs[i], twiddle_d, checksum, bs, thread_bs);
    }

    cudaEventRecord(fft_begin);
    #pragma unroll
    for (int j = 0; j < ntest; ++j){
    
        #pragma unroll
        for(int i = 0; i < kernel_launch_times; ++i){
            entry.turboFFTArr[logN][i]<<<griddims[i], blockdims[i], shared_size[i]>>>(inputs[i], outputs[i], twiddle_d, checksum, bs, thread_bs);
            cudaDeviceSynchronize();
        }
    
        cudaDeviceSynchronize();
    }
    cudaEventRecord(fft_end);
    cudaEventSynchronize(fft_begin);
    cudaEventSynchronize(fft_end);
    cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);
    elapsed_time = elapsed_time / ntest;
    gflops = 5 * N * log2f(N) * bs / elapsed_time * 1000 / 1000000000.f;
    mem_bandwidth = (float)(N * bs * sizeof(DataType) * 2) / (elapsed_time) * 1000.f / 1000000000.f;
    printf("turboFFT, %d, %d, %8.3f, %8.3f, %8.3f\n",  (int)log2f(N),  (int)log2f(bs), elapsed_time, gflops, mem_bandwidth);
    
    checkCudaErrors(cudaMemcpy((void*)output_turbofft, (void*)outputs[kernel_launch_times - 1], N * bs * sizeof(DataType), cudaMemcpyDeviceToHost));
}

template <typename DataType, int if_ft, int if_err>
void TurboFFT_main(ProgramConfig &config){
    DataType* input, *output_turbofft, *output_cufft;
    DataType* input_d, *output_d, *twiddle_d;
    int ntest = 10;

    std::vector<std::vector<long long int>> params;
    
    params = utils::load_parameters(config.param_file_path);

    DataType* checksum_d, *checksum_h;
    cudaMalloc((void**)&checksum_d, sizeof(DataType) * 16384 * 2);
    checksum_h = (DataType*)calloc(16384 * 2, sizeof(DataType));
    DataType* dest = checksum_h;
    for(int i = 2; i <= (1 << 13); i *= 2){
        utils::getDFTMatrixChecksum(dest, i);
        dest += i;
    }
    cudaMemcpy((void*)checksum_d, (void*)checksum_h, sizeof(DataType) * 16384 * 2, cudaMemcpyHostToDevice);
    if(!config.if_bench){
        // Verification
        utils::initializeData<DataType>(input, input_d, output_d, output_turbofft, output_cufft, twiddle_d, config.N, config.bs_end);

        if(config.if_verify){
            test_turbofft<DataType, if_ft, if_err>(input_d, output_d, output_turbofft, twiddle_d, checksum_d, params[config.logN], config.bs, config.thread_bs, 1);
            profiler::cufft::test_cufft<DataType>(input_d, output_d, output_cufft, config.N, config.bs, 1);            
            utils::compareData<DataType>(output_turbofft, output_cufft, config.N * config.bs, 1e-4);
        }
        // Profiling
        if(config.if_profile){
            long long int bs_begin = config.bs;
            for(int bs = bs_begin; bs <= config.bs_end; bs += config.bs_gap)
            profiler::cufft::test_cufft<DataType>(input_d, output_d, output_cufft, config.N, config.bs, ntest);
            
            for(int bs = bs_begin; bs <= config.bs_end; bs += config.bs_gap)
            test_turbofft<DataType, if_ft, if_err>(input_d, output_d, output_turbofft, twiddle_d, checksum_d, params[config.logN], config.bs, config.thread_bs, ntest);
        }
    }
    
    if(config.if_bench){
        utils::initializeData<DataType>(input, input_d, output_d, output_turbofft, output_cufft, twiddle_d, 1 << 25, 16 + 3);
        long long int N = 1;
        for(long long int logN = 1; logN <= 25; ++logN){
            N *= 2;
            long long int bs = 1;
            // bs = bs << (28-logN);
            for(int i = 0; i < 29 - logN; i += 1){
                // profiler::cufft::test_cufft<DataType>(input_d, output_d, output_cufft, N, bs, ntest);
                test_turbofft<DataType, if_ft, if_err>(input_d, output_d, output_turbofft, twiddle_d, checksum_d, params[logN], bs, config.thread_bs, ntest);
                bs *= 2;
                // break; 
            }
        }
    }
    cudaFree(input_d);
    cudaFree(output_d);
    cudaFree(twiddle_d);
    free(input);
    free(output_cufft);
    free(output_turbofft);
}

int main(int argc, char *argv[]){
    ProgramConfig config;
    config.parseCommandLine(argc, argv);
    
    config.displayConfig();
    // Proceed with the rest of the program

    if(config.datatype == 0) {
        if(config.if_ft == 0) TurboFFT_main<float2, 0, 0>(config);
        else if(config.if_err == 0) TurboFFT_main<float2, 1, 0>(config);
        else TurboFFT_main<float2, 1, 1>(config);
    }
    else {
        if(config.if_ft == 0) TurboFFT_main<double2, 0, 0>(config);
        else if(config.if_err == 0) TurboFFT_main<double2, 1, 0>(config);
        else TurboFFT_main<double2, 1, 1>(config);
    }
    
    return 0;
}


    