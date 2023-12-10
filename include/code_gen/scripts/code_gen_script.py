import pandas as pd
def code_gen_script():
    df = pd.read_csv('parameter_radix2_vkfft.csv')
    ft_fft_script = '''
#include <stdlib.h>
#include <complex>
#include "kernels.cuh"
#include <cuda_runtime.h> 
#include <cufftXt.h>
#include "utils/utils.cuh"   
# define RADIX 2
#define FLOAT2_NORM(a, res) res = a.x * a.x + a.y * a.y;
int main(int argc, char** argv){  
    // #if (V == 1)
    int __log_N__, __log_N_st__ = 3;
    float * t_cufft, *t_vkfft, *t_fft;
    t_cufft = (float*)malloc(sizeof(float) * 34);
    t_vkfft = (float*)malloc(sizeof(float) * 34);
    t_fft = (float*)malloc(sizeof(float) * 34);
    
    if (argc < 2){
        printf("Please input log(N)\\n");
        return -1;
    }
    else if(argc == 2) __log_N__ = atoi(argv[1]);
    else if(argc == 3){
        __log_N__ = atoi(argv[2]);
        __log_N_st__ = atoi(argv[1]);
    }
    // #endif
    long long N = pow((double)RADIX, (double)__log_N__); 
    int random_seed = 10;  
    #if P_FFT == 1
    int num_tests = 10;
    #else
    int num_tests = 1;
    #endif
    srandom(random_seed); 
    float *input = (float*)calloc(N * 2, sizeof(float)); 
    float *output_ref, *output;
    
    output_ref = (float*)calloc(N * 2, sizeof(float));
    output = (float*)calloc(N * 2, sizeof(float));
    
    float r[6];
    
    r[0] = 1.0f;
    r[1] = 0.0f;
    r[2] = -0.5f;
    r[3] = -0.8660253882408142f;
    r[4] = -0.5f;
    r[5] = 0.8660253882408142f;
    for(int i = 0; i < 3; ++i){
        r[i * 2] = cosf(-2 * M_PI * (i % 3) / 3);
        r[i * 2 + 1] = sinf(-2 * M_PI * (i % 3) / 3);
    }
    
    float *input_d, *output_d, *output_d_vkfft, *output_d_cufft, *output_d_1, *output_d_ref_1, *checksum_r, *checksum_r_d, *dftmtx;
    checksum_r = (float*)calloc(1024*2, sizeof(float));
    dftmtx = (float*)calloc(1024*1024*2, sizeof(float));
    CUDA_CALLER(cudaMalloc((void**)&input_d, sizeof(float) * N * 2));
    CUDA_CALLER(cudaMalloc((void**)&output_d, sizeof(float) * N * 2));
    
    

    for(int i = 0; i < N * 2; ++i){ 
            input[i] = (float)(random() % 100) / (float)100;
    }
    '''
    r_ = 8
    for i in range(3, 11):
        ft_fft_script += f'''
        float* checksum_r_{i}, *checksum_r_d_{i};
        checksum_r_{i} = (float*)calloc({r_}*2, sizeof(float));
        CUDA_CALLER(cudaMalloc((void**)&checksum_r_d_{i}, sizeof(float) * {r_} * 2));
        // printf("################################### {r_} ###################################\\n");
        for(int i = 0; i < {r_}; ++i)
        '''
        ft_fft_script += '''
        {
        '''
        ft_fft_script += f'''
        for(int j = 0; j < {r_}; ++j )
        '''
        ft_fft_script += '''
        {
        '''
        ft_fft_script += f'''
            dftmtx[i + (j * 2) * {r_}] = cosf((float)(-2 * M_PI * i * j) / {r_}.f);
            dftmtx[i + (j * 2 + 1) * {r_}] = sinf((float)(-2 * M_PI * i * j) / {r_}.f);
        '''
        ft_fft_script += '''
        }
    }
    '''
        ft_fft_script += f'''
    for(int i = 0; i < {r_}; ++i)
    '''
        ft_fft_script += '''
    {
    '''
        ft_fft_script += f'''
        checksum_r_{i}[i * 2] = 0;
        checksum_r_{i}[i * 2 + 1] = 0;
        for(int j = 0; j < {r_}; ++j)
    '''
        ft_fft_script += '''
    {
    '''
        ft_fft_script += f'''
            float real = dftmtx[j + i * 2 * {r_}];
            float imag = dftmtx[j + (i * 2 + 1) * {r_}];
            checksum_r_{i}[i * 2] += real * r[(j % 3) * 2] - imag * r[(j % 3) * 2 + 1];
            checksum_r_{i}[i * 2 + 1] += imag * r[(j % 3) * 2] + real * r[(j % 3) * 2 + 1];
    '''
        ft_fft_script += '''    
    }
    }
    '''
    #     ft_fft_script += f'''
    #     for(int i = 0; i < {r_}; ++i)
    #     printf("%d, %f, %f\\n", i,  checksum_r_{i}[i * 2], checksum_r_{i}[i * 2 + 1] );
    # '''
        ft_fft_script += f'''
    cudaMemcpy((void*)checksum_r_d_{i}, (void*)checksum_r_{i}, 2 * {r_} * sizeof(float), cudaMemcpyHostToDevice);
    '''
        r_ *= 2
    ft_fft_script += '''
    
    cudaMemcpy((void*)input_d, (void*)input, 2 * N * sizeof(float), cudaMemcpyHostToDevice);
    

    cufftHandle plan;  


    cudaEvent_t fft_begin, fft_end;
    float elapsed_time_vkfft, elapsed_time, elapsed_time_cufft; 
    std::chrono::steady_clock::time_point timeSt; // = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point timeEnd; // = std::chrono::steady_clock::now();
	float totTime, totTime_vkfft, totTime_cufft;
    cudaEventCreate(&fft_begin);
    cudaEventCreate(&fft_end);
    
    #if P_FFT == 1
    for(int log_N = __log_N_st__; log_N <= __log_N__; ++log_N){
    #else
    int log_N = __log_N__;
    #endif
    N = pow(double(RADIX), double(log_N));
    '''
    
    N = 3
    
    while N <= 13:
        ft_fft_script += f'''
    if(log_N == {N})''' + '''{
        '''
        for i in range(1,2):
            if  int(df[f'sm_size_{i}'][N-1]) >= 65536:
                ft_fft_script += f'''
        cudaFuncSetAttribute(fft_radix2_logN{int(df['logN'][N-1])}, cudaFuncAttributeMaxDynamicSharedMemorySize, {int(df[f'sm_size_{i}'][N-1])});
        // cudaFuncSetAttribute(VkFFT_main_logN{int(df['logN'][N-1])}, cudaFuncAttributeMaxDynamicSharedMemorySize, {int(df[f'sm_size_{i}'][N-1])});
        '''
        ft_fft_script += '''
        cudaEventCreate(&fft_begin);
        cudaEventCreate(&fft_end);
        '''
        # ft_fft_script += '''
        # {
        # '''
        # ft_fft_script += f'''
        #     dim3 gridDim({int(df['num_block_1'][N-1])}, 1, 1);
        #     dim3 blockDim({int(df['blockdim_x_1'][N-1])}, {int(df['blockdim_y_1'][N-1])}, 1);
        #     cudaEventRecord(fft_begin);
        #     timeSt = std::chrono::steady_clock::now();
        # '''
        # ft_fft_script += '''
        #     for(int i = 0; i < num_tests; ++i){
        # '''
        # ft_fft_script += f'''
        #         VkFFT_main_logN{int(df['logN'][N-1])} <<<gridDim, blockDim, {int(df['sm_size_1'][N-1])}>>>((float2*)input_d, (float2*)output_d);
        #         cudaDeviceSynchronize();  
        # '''
        # ft_fft_script += '''
        #     }
        #     timeEnd = std::chrono::steady_clock::now();
        #     totTime_vkfft = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
        #     cudaEventRecord(fft_end);
        #     cudaEventSynchronize(fft_begin);  
        #     cudaEventSynchronize(fft_end);
        #     cudaEventElapsedTime(&elapsed_time_vkfft, fft_begin, fft_end);
        # }
        ft_fft_script += '''
        {
            cufftCreate(&plan);
            cufftPlan1d(&plan, N, CUFFT_C2C, 1); 
            cudaEventRecord(fft_begin);
            timeSt = std::chrono::steady_clock::now();
            for(int i = 0; i < num_tests; ++i){
                cufftExecC2C(plan, (cufftComplex *)input_d, (cufftComplex *)output_d, CUFFT_FORWARD);
                cudaDeviceSynchronize(); 
            } 
            timeEnd = std::chrono::steady_clock::now();
            totTime_cufft = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
            cudaEventRecord(fft_end);  
            cudaEventSynchronize(fft_begin);
            cudaEventSynchronize(fft_end);
            cudaEventElapsedTime(&elapsed_time_cufft, fft_begin, fft_end);  
            cudaMemcpy((void*)output_ref, (void*)output_d, 2 * N * sizeof(float), cudaMemcpyDeviceToHost); 
            cufftDestroy(plan);
        }
    '''
        ft_fft_script += '''
        {
        '''
        ft_fft_script += f'''
            CUDA_CALLER(cudaMalloc((void**)&output_d_1, sizeof(float) * N * 2));
            dim3 gridDim({int(df['num_block_1'][N-1])}, 1, 1);
            dim3 blockDim({int(df['blockdim_x_1'][N-1])}, {int(df['blockdim_y_1'][N-1])}, 1);
            cudaEventRecord(fft_begin);
            timeSt = std::chrono::steady_clock::now();
            '''
        ft_fft_script += '''
            for(int i = 0; i < num_tests; ++i){
        '''
        ft_fft_script += f'''
                fft_radix2_logN{int(df['logN'][N-1])} <<<gridDim, blockDim, {int(df['sm_size_1'][N-1])}>>> ((float2*)input_d, (float2*)output_d);
                cudaDeviceSynchronize();
        '''
        ft_fft_script += '''
            }
            timeEnd = std::chrono::steady_clock::now();
            totTime = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
            cudaEventRecord(fft_end);  
            cudaEventSynchronize(fft_begin);
            cudaEventSynchronize(fft_end);
            cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);
            cudaMemcpy((void*)output, (void*)output_d, 2 * N * sizeof(float), cudaMemcpyDeviceToHost);
            CUDA_CALLER(cudaFree(output_d_1));
        }
    }
        '''
        N += 1
    
    # N += 23
    while N <= 22:
        ft_fft_script += f'''
    if(log_N == {N})''' + '''{
        '''
        for i in range(1,3):
            if  int(df[f'sm_size_{i}'][N-1]) >= 65536:
                ft_fft_script += f'''
        cudaFuncSetAttribute(fft_radix2_logN{int(df['logN'][N-1])}_{i}, cudaFuncAttributeMaxDynamicSharedMemorySize, {int(df[f'sm_size_{i}'][N-1])});
        // cudaFuncSetAttribute(VkFFT_main_logN{int(df['logN'][N-1])}_{i}, cudaFuncAttributeMaxDynamicSharedMemorySize, {int(df[f'sm_size_{i}'][N-1])});
        '''
        ft_fft_script += '''
        cudaEventCreate(&fft_begin);
        cudaEventCreate(&fft_end);
        {
        '''
        ft_fft_script += f'''
            CUDA_CALLER(cudaMalloc((void**)&output_d_1, sizeof(float) * N * 2));
            cudaEventRecord(fft_begin);
            timeSt = std::chrono::steady_clock::now();
            '''
        ft_fft_script += '''
            for(int i = 0; i < num_tests; ++i){
        '''
        ft_fft_script += f'''{{
                dim3 gridDim({int(df['num_block_1'][N-1])}, 1, 1);
                dim3 blockDim({int(df['blockdim_x_1'][N-1])}, {int(df['blockdim_y_1'][N-1])}, 1);
                fft_radix2_logN{int(df['logN'][N-1])}_1 <<<gridDim, blockDim, {int(df['sm_size_1'][N-1])}>>> ((float2*)input_d, (float2*)output_d_1);
                cudaDeviceSynchronize();
            }}
        '''
        ft_fft_script += f'''{{
                dim3 gridDim({int(df['num_block_2'][N-1])}, 1, 1);
                dim3 blockDim({int(df['blockdim_x_2'][N-1])}, {int(df['blockdim_y_2'][N-1])}, 1);
                fft_radix2_logN{int(df['logN'][N-1])}_2 <<<gridDim, blockDim, {int(df['sm_size_2'][N-1])}>>> ((float2*)output_d_1, (float2*)output_d);
                cudaDeviceSynchronize();
            }}
        '''
        ft_fft_script += '''
            }
            timeEnd = std::chrono::steady_clock::now();
            totTime = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
            cudaEventRecord(fft_end);  
            cudaEventSynchronize(fft_begin);
            cudaEventSynchronize(fft_end);
            cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);
            cudaMemcpy((void*)output, (void*)output_d, 2 * N * sizeof(float), cudaMemcpyDeviceToHost);
            CUDA_CALLER(cudaFree(output_d_1));
        }
        '''
        # ft_fft_script += '''
        # {
        # '''
        # ft_fft_script += f'''
        #     cudaEventRecord(fft_begin);
        #     timeSt = std::chrono::steady_clock::now();
        # '''
        # ft_fft_script += '''
        #     for(int i = 0; i < num_tests; ++i){
        # '''
        # ft_fft_script += f'''{{
        #         dim3 gridDim({int(df['num_block_1'][N-1])}, 1, 1);
        #         dim3 blockDim({int(df['blockdim_x_1'][N-1])}, {int(df['blockdim_y_1'][N-1])}, 1);
        #         VkFFT_main_logN{int(df['logN'][N-1])}_1 <<<gridDim, blockDim, {int(df['sm_size_1'][N-1])}>>>((float2*)input_d, (float2*)output_d_1);
        #         cudaDeviceSynchronize();  
        # }}
        # '''
        # ft_fft_script += f'''{{
        #         dim3 gridDim({int(df['num_block_2'][N-1])}, 1, 1);
        #         dim3 blockDim({int(df['blockdim_x_2'][N-1])}, {int(df['blockdim_y_2'][N-1])}, 1);
        #         VkFFT_main_logN{int(df['logN'][N-1])}_2 <<<gridDim, blockDim, {int(df['sm_size_2'][N-1])}>>>((float2*)output_d_1, (float2*)output_d);
        #         cudaDeviceSynchronize();  
        # }}
        # '''
        # ft_fft_script += '''
        #     }
        #     timeEnd = std::chrono::steady_clock::now();
        #     totTime_vkfft = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
        #     cudaEventRecord(fft_end);
        #     cudaEventSynchronize(fft_begin);  
        #     cudaEventSynchronize(fft_end);
        #     cudaEventElapsedTime(&elapsed_time_vkfft, fft_begin, fft_end);
        # }
        ft_fft_script += '''
        {
            cufftCreate(&plan);
            cufftPlan1d(&plan, N, CUFFT_C2C, 1); 
            cudaEventRecord(fft_begin);
            timeSt = std::chrono::steady_clock::now();
            for(int i = 0; i < num_tests; ++i){
                cufftExecC2C(plan, (cufftComplex *)input_d, (cufftComplex *)output_d, CUFFT_FORWARD);
                cudaDeviceSynchronize(); 
            } 
            timeEnd = std::chrono::steady_clock::now();
            totTime_cufft = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
            cudaEventRecord(fft_end);  
            cudaEventSynchronize(fft_begin);
            cudaEventSynchronize(fft_end);
            cudaEventElapsedTime(&elapsed_time_cufft, fft_begin, fft_end);   
            cudaMemcpy((void*)output_ref, (void*)output_d, 2 * N * sizeof(float), cudaMemcpyDeviceToHost);
            cufftDestroy(plan);
        }
    }    
    '''
        N += 1
    
    while N <= 29:
        ft_fft_script += f'''
    if(log_N == {N})''' + '''{
        '''
        for i in range(1,4):
            if  int(df[f'sm_size_{i}'][N-1]) >= 65536:
                ft_fft_script += f'''
        cudaFuncSetAttribute(fft_radix2_logN{int(df['logN'][N-1])}_{i}, cudaFuncAttributeMaxDynamicSharedMemorySize, {int(df[f'sm_size_{i}'][N-1])});
        // cudaFuncSetAttribute(VkFFT_main_logN{int(df['logN'][N-1])}_{i}, cudaFuncAttributeMaxDynamicSharedMemorySize, {int(df[f'sm_size_{i}'][N-1])});
        '''
        # ft_fft_script += '''
        # {
        # cudaEventCreate(&fft_begin);
        # cudaEventCreate(&fft_end);
        # '''
        # ft_fft_script += f'''
        #     cudaEventRecord(fft_begin);
        #     timeSt = std::chrono::steady_clock::now();
        # '''
        # ft_fft_script += '''
        #     for(int i = 0; i < num_tests; ++i){
        # '''
        # ft_fft_script += f'''{{
        #         dim3 gridDim({int(df['num_block_1'][N-1])}, 1, 1);
        #         dim3 blockDim({int(df['blockdim_x_1'][N-1])}, {int(df['blockdim_y_1'][N-1])}, 1);
        #         VkFFT_main_logN{int(df['logN'][N-1])}_1 <<<gridDim, blockDim, {int(df['sm_size_1'][N-1])}>>>((float2*)input_d, (float2*)output_d);
        #         cudaDeviceSynchronize();  
        # }}
        # '''
        # ft_fft_script += f'''{{
        #         dim3 gridDim({int(df['num_block_2'][N-1])}, 1, 1);
        #         dim3 blockDim({int(df['blockdim_x_2'][N-1])}, {int(df['blockdim_y_2'][N-1])}, 1);
        #         VkFFT_main_logN{int(df['logN'][N-1])}_2 <<<gridDim, blockDim, {int(df['sm_size_2'][N-1])}>>>((float2*)output_d, (float2*)output_d_1);
        #         cudaDeviceSynchronize();  
        # }}
        # '''
        # ft_fft_script += f'''{{
        #         dim3 gridDim({int(df['num_block_3'][N-1])}, 1, 1);
        #         dim3 blockDim({int(df['blockdim_x_3'][N-1])}, {int(df['blockdim_y_3'][N-1])}, 1);
        #         VkFFT_main_logN{int(df['logN'][N-1])}_3 <<<gridDim, blockDim, {int(df['sm_size_3'][N-1])}>>>((float2*)output_d_1, (float2*)output_d);
        #         cudaDeviceSynchronize();  
        # }}
        # '''
        # ft_fft_script += '''
        #     }
        #     timeEnd = std::chrono::steady_clock::now();
        #     totTime_vkfft = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
        #     cudaEventRecord(fft_end);
        #     cudaEventSynchronize(fft_begin);  
        #     cudaEventSynchronize(fft_end);
        #     cudaEventElapsedTime(&elapsed_time_vkfft, fft_begin, fft_end);
        # }
        ft_fft_script += '''
        {
            cufftCreate(&plan);
            cudaEventCreate(&fft_begin);
            cudaEventCreate(&fft_end);
            int res = cufftPlan1d(&plan, N, CUFFT_C2C, 1); 
            if(res != 0){
                printf("cuFFT plan error!\\n");
                return 0; 
            } 
            cudaEventRecord(fft_begin);
            timeSt = std::chrono::steady_clock::now();
            for(int i = 0; i < num_tests; ++i){
                cufftExecC2C(plan, (cufftComplex *)input_d, (cufftComplex *)output_d, CUFFT_FORWARD);
                cudaDeviceSynchronize(); 
            } 
            timeEnd = std::chrono::steady_clock::now();
            totTime_cufft = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
            cudaEventRecord(fft_end);  
            cudaEventSynchronize(fft_begin);
            cudaEventSynchronize(fft_end);
            cudaEventElapsedTime(&elapsed_time_cufft, fft_begin, fft_end);   
            
            CUDA_CALLER(cudaMemcpy((void*)output_ref, (void*)output_d, 2 * N * sizeof(float), cudaMemcpyDeviceToHost));
            cufftDestroy(plan);

        }    
        
    '''
        ft_fft_script += '''
        {
        cudaEventCreate(&fft_begin);
        cudaEventCreate(&fft_end);
        float * reduction = (float*)calloc(2 * 65536, sizeof(float));
        float * reduction_d, *global_checksum_d;
        CUDA_CALLER(cudaMalloc((void**)&reduction_d, sizeof(float) * 65536 * 2));;
        CUDA_CALLER(cudaMalloc((void**)&global_checksum_d, sizeof(float) * 2));;
        cudaMemcpy((void*)reduction_d, (void*)reduction, 2 * 65536 * sizeof(float), cudaMemcpyHostToDevice);
        CUDA_CALLER(cudaMalloc((void**)&output_d_1, sizeof(float) * N * 2));
        cublasHandle_t handle;
        cublasCreate(&handle);
        '''
        
        ft_fft_script += f'''
            cudaEventRecord(fft_begin);
            timeSt = std::chrono::steady_clock::now();
            '''
        ft_fft_script += '''
            for(int i = 0; i < num_tests; ++i){
        '''
        ft_fft_script += f'''{{
                dim3 gridDim({int(df['num_block_1'][N-1])}, 1, 1);
                dim3 blockDim({int(df['blockdim_x_1'][N-1])}, {int(df['blockdim_y_1'][N-1])}, 1);
                fft_radix2_logN{int(df['logN'][N-1])}_1 <<<gridDim, blockDim, {int(df['sm_size_1'][N-1])}>>> ((float2*)input_d, (float2*)output_d, (float2*) checksum_r_d_{int(df['logN1'][N-1])});
                cudaDeviceSynchronize();
            }}
            #if defined(GLOBAL_ON)
            cublasSdot(handle, {int(df['num_block_3'][N-1])}, reduction_d, 1, reduction_d + {int(df['num_block_3'][N-1])}, 1, global_checksum_d);
            #endif
        '''
        ft_fft_script += f'''{{
                dim3 gridDim({int(df['num_block_2'][N-1])}, 1, 1);
                dim3 blockDim({int(df['blockdim_x_2'][N-1])}, {int(df['blockdim_y_2'][N-1])}, 1);
                fft_radix2_logN{int(df['logN'][N-1])}_2 <<<gridDim, blockDim, {int(df['sm_size_2'][N-1])}>>> ((float2*)output_d, (float2*)output_d_1, (float2*) checksum_r_d_{int(df['logN2'][N-1])});
                cudaDeviceSynchronize();
            }}
            #if defined(GLOBAL_ON)
            cublasSdot(handle, {int(df['num_block_3'][N-1])}, reduction_d, 1, reduction_d + {int(df['num_block_3'][N-1])}, 1, global_checksum_d);
            #endif
        '''
            
        ft_fft_script += f'''{{
                dim3 gridDim({int(df['num_block_3'][N-1])}, 1, 1);
                dim3 blockDim({int(df['blockdim_x_3'][N-1])}, {int(df['blockdim_y_3'][N-1])}, 1);
                fft_radix2_logN{int(df['logN'][N-1])}_3 <<<gridDim, blockDim, {int(df['sm_size_3'][N-1])}>>> ((float2*)output_d_1, (float2*)output_d, (float2*) checksum_r_d_{int(df['logN3'][N-1])});
                cudaDeviceSynchronize();
            }}
        '''
        ft_fft_script += '''
            
            }
            timeEnd = std::chrono::steady_clock::now();
            totTime = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeSt).count();
            cudaEventRecord(fft_end);  
            cudaEventSynchronize(fft_begin);
            cudaEventSynchronize(fft_end);
            
            cudaEventElapsedTime(&elapsed_time, fft_begin, fft_end);
            cudaFree(reduction_d);
            cudaFree(global_checksum_d);
            CUDA_CALLER(cudaMemcpy((void*)output, (void*)output_d, 2 * N * sizeof(float), cudaMemcpyDeviceToHost));
            CUDA_CALLER(cudaFree(output_d_1));
        }
        }
        '''
        N += 1
    
    ft_fft_script += '''
    #if V_FFT == 1
    // cudaMemcpy((void*)output_ref, (void*)output_d_vkfft, 2 * N * sizeof(float), cudaMemcpyDeviceToHost);
    // cudaMemcpy((void*)output, (void*)output_d, sizeof(float) * 2 * N, cudaMemcpyDeviceToHost);
    // cudaMemcpy((void*)output, (void*)output_d_cufft, sizeof(float) * 2 * N, cudaMemcpyDeviceToHost);
    // cudaMemcpy((void*)output, (void*)output_d_ref_1, 2 * N * sizeof(float), cudaMemcpyDeviceToHost);
    // cudaMemcpy((void*)output_ref, (void*)output_d_1, sizeof(float) * 2 * N, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    bool pass = true;
    for(int i = 0; i < 2 * N; i +=2){
        float2 res = *(float2*)(output + i); 
        float2 res_ref = *(float2*)(output_ref + i);
        float norm, norm_ref; 
        FLOAT2_NORM(res, norm);
        FLOAT2_NORM(res_ref, norm_ref);
        
        float err = fabs(norm - norm_ref);
        if(i % 1000000 ==0){
        printf("error %f detected at %d\\n", err / fabs(norm), i / 2);
        printf("ref[%d]: %.3f + %.3f i\\n",  i / 2, res_ref.x, res_ref.y);
        printf("res[%d]: %.3f + %.3f i\\n\\n",  i / 2, res.x, res.y);
        }
        if(err / fabs(norm) > 0.05){
            printf("error %f detected at %d\\n", err / fabs(norm), i / 2);
            printf("ref[%d]: %.3f + %.3f i\\n",  i / 2, res_ref.x, res_ref.y);
            printf("res[%d]: %.3f + %.3f i\\n\\n",  i / 2, res.x, res.y);
            pass = false;
            // break;
        }   
    }
    if(pass) printf("Pass!\\n");
    else printf("Fail!\\n");
    #endif

    #if P_FFT == 1
    elapsed_time /= num_tests;
    elapsed_time_vkfft /= num_tests;
    elapsed_time_cufft /= num_tests;
    totTime /= num_tests;
    totTime_vkfft /= num_tests;
    totTime_cufft /= num_tests;
    if(log_N == __log_N_st__)printf("| SIZE |  Execution Time (us)             |   Shared   | #threads |\\n");
    if(log_N == __log_N_st__)printf("|log(N)|   Ours   |   cuFFT   | Memory (KB)|          |\\n");
    // if(log_N == __log_N_st__)printf("|log(N)|   Ours   |   VkFFT   |   cuFFT   | Memory (KB)|          |\\n");
    // printf("|%6d| %8.3f | %8.3f  |%8.3f   |%8.3f    |%10d|\\n", int(log2f((float)N)), elapsed_time * 1000, elapsed_time_vkfft * 1000, elapsed_time_cufft * 1000, (float)sizeof(float) * (float)N * 2.f / 1024.f, N / 8);
    printf("|%6d| %8.3f | %8.3f  |%8.3f   |%8.3f    |%10d|\\n", int(log2f((float)N)), elapsed_time * 1000, elapsed_time_cufft * 1000, (float)sizeof(float) * (float)N * 2.f / 1024.f, N / 8);
    t_fft[log_N] = elapsed_time;
    t_cufft[log_N] = elapsed_time_cufft;
    t_vkfft[log_N] = elapsed_time_vkfft;
    }
    printf("Execution Time\\n");
    printf("t_fft = th.as_tensor([");
    for(int i = __log_N_st__; i <= __log_N__; ++i ){
        printf("%8f,", t_fft[i]);
    }
    printf("])\\n");

    printf("t_cufft = th.as_tensor([");
    for(int i = __log_N_st__; i <= __log_N__; ++i ){
        printf("%8f,", t_cufft[i]);
    }
    printf("])\\n");
    
    // printf("t_vkfft = th.as_tensor([");
    // for(int i = __log_N_st__; i <= __log_N__; ++i ){
    //     printf("%8f,", t_vkfft[i]);
    // }
    // printf("])\\n");

    printf("\\n Flops\\n");
    printf("gflops_fft = th.as_tensor([");
    for(int i = __log_N_st__; i <= __log_N__; ++i ){
        long long N = pow((double)RADIX, (double)i);
        printf("%8.1f,", 5 * N * i / t_fft[i] * 1000.f / 1000000000.f);
    }
    printf("])\\n");

    printf("gflops_cufft = th.as_tensor([");
    for(int i = __log_N_st__; i <= __log_N__; ++i ){
        long long N = pow((double)RADIX, (double)i);
        printf("%8.1f,", 5 * N * i / t_cufft[i] * 1000.f / 1000000000.f);
    }
    printf("])\\n");
    
    // printf("gflops_vkfft = th.as_tensor([");
    // for(int i = __log_N_st__; i <= __log_N__; ++i ){
    //     long long N = pow((double)RADIX, (double)i);
    //     printf("%8.1f,", 5 * N * i / t_vkfft[i] * 1000.f / 1000000000.f);
    // }
    // printf("])\\n");
    #endif
    return 0;
}

    '''
    return ft_fft_script