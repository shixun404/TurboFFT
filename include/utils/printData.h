#include "utils.h"
namespace utils{
template <typename DataType>
void printData(DataType* res, long long int N){
    double rel_error = 0.;
    for(int i = 0; i < N; ++i){
        printf("res[%d] = %f + %f j\n", i, res[i].x, res[i].y);
    }
}

template <>
void printData<nv_bfloat162>(nv_bfloat162* input, long long int N){

    double rel_error = 0.;
    for(int i = 0; i < N; ++i){
        float2 input_i =  __bfloat1622float2(input[i]);
        printf("res[%d] = %f + %f j\n", i, input_i.x, input_i.y);
    }
}

template <>
void printData<half2>(half2* input, long long int N){

    double rel_error = 0.;
    for(int i = 0; i < N; ++i){
        float2 input_i =  __half22float2(input[i]);
        printf("res[%d] = %f + %f j\n", i, input_i.x, input_i.y);
    }
}
}