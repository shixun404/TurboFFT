#include "utils.h"
namespace utils{
    template <typename DataType>
    void initializeData(DataType *&input, DataType *&input_d, DataType *&output_d, 
                    DataType *&output_turbofft, DataType *&output_cufft, long long int N){
    input = (DataType*)calloc(N, sizeof(DataType));
    output_turbofft = (DataType*)calloc(N, sizeof(DataType));
    output_cufft = (DataType*)calloc(N, sizeof(DataType));
    checkCudaErrors(cudaMalloc((void**)&input_d, sizeof(DataType) * N));
    checkCudaErrors(cudaMalloc((void**)&output_d, sizeof(DataType) * N));
    
    for(int i = 0; i < N; ++i){
        input[i].x = 1;
        input[i].y = 1;
    }

    checkCudaErrors(cudaMemcpy((void*)input_d, (void*)input, N * sizeof(DataType), cudaMemcpyHostToDevice));
    }

    template <>
    void initializeData<nv_bfloat162>(nv_bfloat162 *&input, nv_bfloat162 *&input_d, 
                                    nv_bfloat162 *&output_d, nv_bfloat162 *&output_turbofft,
                                    nv_bfloat162 *&output_cufft, long long int N)
    {
        input = (nv_bfloat162*)calloc(N, sizeof(nv_bfloat162));
        output_turbofft = (nv_bfloat162*)calloc(N, sizeof(nv_bfloat162));
        output_cufft = (nv_bfloat162*)calloc(N, sizeof(nv_bfloat162));
        checkCudaErrors(cudaMalloc((void**)&input_d, sizeof(nv_bfloat162) * N));
        checkCudaErrors(cudaMalloc((void**)&output_d, sizeof(nv_bfloat162) * N));
        
        for(int i = 0; i < N; ++i){
            float2 input_i_f32 = {1.0f, 1.0f};
            input[i] = __float22bfloat162_rn(input_i_f32);
        }

        checkCudaErrors(cudaMemcpy((void*)input_d, (void*)input, N * sizeof(nv_bfloat162), cudaMemcpyHostToDevice));
    }
}
