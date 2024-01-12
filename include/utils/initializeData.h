#include "utils.h"
namespace utils{
    template <typename DataType>
    void initializeData(DataType *&input, DataType *&input_d, DataType *&output_d, 
                    DataType *&output_turbofft, DataType *&output_cufft, long long int N, long long int bs){
    bs = bs;
    input = (DataType*)calloc(N * bs, sizeof(DataType));
    output_turbofft = (DataType*)calloc(N * bs, sizeof(DataType));
    output_cufft = (DataType*)calloc(N * bs, sizeof(DataType));
    int res = cudaMalloc((void**)&input_d, sizeof(DataType) * N * bs);
    
    printf("%lld, %lld, %lld, %lld\n", sizeof(DataType), N, bs,
    (long long int)sizeof(DataType) * (long long int)N * (long long int)bs);
    printf("Intiliaze input_d status %d\n", res);
    if(res) exit(-1);
    // checkCudaErrors(cudaMalloc((void**)&output_d, sizeof(DataType) * N * bs));
    res = cudaMalloc((void**)&output_d, sizeof(DataType) * N * bs);
    printf("Intiliaze output_d status %d\n", res);
    if(res) exit(-1);
    for(int i = 0; i < N * bs; ++i){
        input[i].x = 1;
        input[i].y = 1;
    }

    checkCudaErrors(cudaMemcpy((void*)input_d, (void*)input, N * bs * sizeof(DataType), cudaMemcpyHostToDevice));
    }

    template <>
    void initializeData<nv_bfloat162>(nv_bfloat162 *&input, nv_bfloat162 *&input_d, 
                                    nv_bfloat162 *&output_d, nv_bfloat162 *&output_turbofft,
                                    nv_bfloat162 *&output_cufft, long long int N, long long int bs)
    {
        input = (nv_bfloat162*)calloc(N, sizeof(nv_bfloat162));
        output_turbofft = (nv_bfloat162*)calloc(N * bs, sizeof(nv_bfloat162));
        output_cufft = (nv_bfloat162*)calloc(N * bs, sizeof(nv_bfloat162));
        checkCudaErrors(cudaMalloc((void**)&input_d, sizeof(nv_bfloat162) * N * bs));
        checkCudaErrors(cudaMalloc((void**)&output_d, sizeof(nv_bfloat162) * N * bs));
        
        for(int i = 0; i < N * bs; ++i){
            float2 input_i_f32 = {1.0f, 1.0f};
            input[i] = __float22bfloat162_rn(input_i_f32);
        }

        checkCudaErrors(cudaMemcpy((void*)input_d, (void*)input, N * bs * sizeof(nv_bfloat162), cudaMemcpyHostToDevice));
    }
}
