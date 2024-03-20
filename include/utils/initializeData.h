#include "utils.h"
namespace utils{
        template <typename DataType>
    void initializeData(DataType *&input, DataType *&input_d, DataType *&output_d, 
                    DataType *&output_turbofft, DataType *&output_cufft, DataType *&twiddle_d, long long int N, long long int bs){
    int random_seed = 10;  
    srandom(random_seed); 
    bs = bs;
    input = (DataType*)calloc(N * bs, sizeof(DataType));
    DataType* twiddle = (DataType*)calloc(N * 2, sizeof(DataType));
    output_turbofft = (DataType*)calloc(N * bs, sizeof(DataType));
    output_cufft = (DataType*)calloc(N * bs, sizeof(DataType));
    int res = cudaMalloc((void**)&input_d, sizeof(DataType) * N * bs);
    res = cudaMalloc((void**)&twiddle_d, sizeof(DataType) * N * 2);
    
    // printf("%lld, %lld, %lld, %lld\n", sizeof(DataType), N, bs,
    // (long long int)sizeof(DataType) * (long long int)N * (long long int)bs);
    // printf("Intiliaze input_d status %d\n", res);
    if(res) exit(-1);
    // checkCudaErrors(cudaMalloc((void**)&output_d, sizeof(DataType) * N * bs));
    res = cudaMalloc((void**)&output_d, sizeof(DataType) * N * bs * 2);
    // printf("Intiliaze output_d status %d\n", res);
    if(res) exit(-1);
    for(long long int i = 0; i < N * bs; ++i){
        input[i].x = (double)(random() % 100 - 50) / (double)1000;
        input[i].y = (double)(random() % 100 - 50) / (double)1000;
        // input[i].x = 1;
        // input[i].y = 0;
    }

    long long int cur_N = 1;
    for(long long int i = 0; i < N * 2 - 1 ; ++i){   
        if(cur_N * 2 - 1 == i) cur_N *= 2;
        double angle = -2.0 * M_PI * (float)(i - cur_N + 1) / (float)cur_N;
        twiddle[i].x = cos(angle);
        twiddle[i].y = sin(angle);
    }

    checkCudaErrors(cudaMemcpy((void*)input_d, (void*)input, N * bs * sizeof(DataType), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy((void*)twiddle_d, (void*)twiddle, N * 2 * sizeof(DataType), cudaMemcpyHostToDevice));
    }

    template <>
    void initializeData<float2>(float2 *&input, float2 *&input_d, float2 *&output_d, 
                    float2 *&output_turbofft, float2 *&output_cufft, float2 *&twiddle_d, long long int N, long long int bs){
    int random_seed = 10;  
    srandom(random_seed); 
    bs = bs;
    input = (float2*)calloc(N * bs, sizeof(float2));
    float2* twiddle = (float2*)calloc(N * 2, sizeof(float2));
    output_turbofft = (float2*)calloc(N * bs, sizeof(float2));
    output_cufft = (float2*)calloc(N * bs, sizeof(float2));
    int res = cudaMalloc((void**)&input_d, sizeof(float2) * N * bs);
    res = cudaMalloc((void**)&twiddle_d, sizeof(float2) * N * 2);
    
    // printf("%lld, %lld, %lld, %lld\n", sizeof(float2), N, bs,
    // (long long int)sizeof(float2) * (long long int)N * (long long int)bs);
    // printf("Intiliaze input_d status %d\n", res);
    if(res) exit(-1);
    // checkCudaErrors(cudaMalloc((void**)&output_d, sizeof(float2) * N * bs));
    res = cudaMalloc((void**)&output_d, sizeof(float2) * N * bs * 2);
    // printf("Intiliaze output_d status %d\n", res);
    if(res) exit(-1);
    for(long long int i = 0; i < N * bs; ++i){
        input[i].x = (double)(random() % 100 - 50) / (double)100;
        input[i].y = (double)(random() % 100 - 50) / (double)100;
    }

    long long int cur_N = 1;
    for(long long int i = 0; i < N * 2 - 1 ; ++i){   
        if(cur_N * 2 - 1 == i) cur_N *= 2;
        double angle = -2.0 * M_PI * (float)(i - cur_N + 1) / (float)cur_N;
        twiddle[i].x = cos(angle);
        twiddle[i].y = sin(angle);
    }

    checkCudaErrors(cudaMemcpy((void*)input_d, (void*)input, N * bs * sizeof(float2), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy((void*)twiddle_d, (void*)twiddle, N * 2 * sizeof(float2), cudaMemcpyHostToDevice));
    }

    // template <>
    // void initializeData<nv_bfloat162>(nv_bfloat162 *&input, nv_bfloat162 *&input_d, 
    //                                 nv_bfloat162 *&output_d, nv_bfloat162 *&output_turbofft,
    //                                 nv_bfloat162 *&output_cufft, long long int N, long long int bs)
    // {
    //     input = (nv_bfloat162*)calloc(N, sizeof(nv_bfloat162));
    //     output_turbofft = (nv_bfloat162*)calloc(N * bs, sizeof(nv_bfloat162));
    //     output_cufft = (nv_bfloat162*)calloc(N * bs, sizeof(nv_bfloat162));
    //     checkCudaErrors(cudaMalloc((void**)&input_d, sizeof(nv_bfloat162) * N * bs));
    //     checkCudaErrors(cudaMalloc((void**)&output_d, sizeof(nv_bfloat162) * N * bs));
        
    //     for(int i = 0; i < N * bs; ++i){
    //         float2 input_i_f32 = {1.0f, 1.0f};
    //         input[i] = __float22bfloat162_rn(input_i_f32);
    //     }

    //     checkCudaErrors(cudaMemcpy((void*)input_d, (void*)input, N * bs * sizeof(nv_bfloat162), cudaMemcpyHostToDevice));
    // }
}
