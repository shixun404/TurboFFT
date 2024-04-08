#include <stdio.h>
#include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cuda_fp16.h>
#include <cuda_bf16.h>
#include <fstream>

#include "utils/utils.h"
#include "turbofft/constants.h"
#include "turbofft/macro_ops.h"
#include "cufft/cufft.h"
#include "cufft/cufft_ft.h"
#include "utils/compareData.h"
#include "utils/printData.h"
#include "utils/initializeData.h"
#include "utils/readCSV.h"
#include "utils/abft.h"
#include "utils/CommandLineParser.h"



template <typename DataType, int if_ft, int if_err_injection>
struct TurboFFT_Kernel_Entry {
void (*turboFFTArr[26][3])(DataType *, DataType *, DataType *, DataType*, int, int);
};

template <typename DataType>
void test_turbofft( DataType* input_d, DataType* output_d, DataType* output_turbofft,
                    DataType* twiddle_d, DataType* checksum, std::vector<long long int> param, 
                    long long int bs, int ntest, auto alpha);



#include "code_gen/generated/float2/fft_radix_2_logN_1_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_2_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_3_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_4_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_5_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_6_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_7_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_8_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_9_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_10_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_11_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_12_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_13_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_13_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_14_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_14_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_15_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_15_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_16_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_16_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_17_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_17_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_18_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_18_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_19_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_19_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_20_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_20_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_21_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_21_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_22_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_22_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_23_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_23_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_23_upload_2.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_24_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_24_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_24_upload_2.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_25_upload_0.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_25_upload_1.cuh"
#include "code_gen/generated/float2/fft_radix_2_logN_25_upload_2.cuh"

#include "code_gen/generated/double2/fft_radix_2_logN_1_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_2_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_3_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_4_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_5_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_6_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_7_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_8_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_9_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_10_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_11_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_12_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_13_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_14_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_14_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_15_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_15_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_16_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_16_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_17_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_17_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_18_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_18_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_19_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_19_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_20_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_20_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_21_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_21_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_22_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_22_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_23_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_23_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_23_upload_2.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_24_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_24_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_24_upload_2.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_25_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_25_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_25_upload_2.cuh"
template<> struct TurboFFT_Kernel_Entry<float2, 0, 0>
{
void (*turboFFTArr [26][3])(float2 *, float2 *, float2 *, float2 *, int, int) ={
 {NULL, NULL, NULL},
{fft_radix_2<float2, 1, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 2, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 3, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 4, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 5, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 6, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 7, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 8, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 9, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 10, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 11, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 12, 0, 0, 0>, NULL, NULL},
{fft_radix_2<float2, 13, 0, 0, 0>, fft_radix_2<float2, 13, 1, 0, 0>, NULL},
{fft_radix_2<float2, 14, 0, 0, 0>, fft_radix_2<float2, 14, 1, 0, 0>, NULL},
{fft_radix_2<float2, 15, 0, 0, 0>, fft_radix_2<float2, 15, 1, 0, 0>, NULL},
{fft_radix_2<float2, 16, 0, 0, 0>, fft_radix_2<float2, 16, 1, 0, 0>, NULL},
{fft_radix_2<float2, 17, 0, 0, 0>, fft_radix_2<float2, 17, 1, 0, 0>, NULL},
{fft_radix_2<float2, 18, 0, 0, 0>, fft_radix_2<float2, 18, 1, 0, 0>, NULL},
{fft_radix_2<float2, 19, 0, 0, 0>, fft_radix_2<float2, 19, 1, 0, 0>, NULL},
{fft_radix_2<float2, 20, 0, 0, 0>, fft_radix_2<float2, 20, 1, 0, 0>, NULL},
{fft_radix_2<float2, 21, 0, 0, 0>, fft_radix_2<float2, 21, 1, 0, 0>, NULL},
{fft_radix_2<float2, 22, 0, 0, 0>, fft_radix_2<float2, 22, 1, 0, 0>, NULL},
{fft_radix_2<float2, 23, 0, 0, 0>, fft_radix_2<float2, 23, 1, 0, 0>, fft_radix_2<float2, 23, 2, 0, 0>},
{fft_radix_2<float2, 24, 0, 0, 0>, fft_radix_2<float2, 24, 1, 0, 0>, fft_radix_2<float2, 24, 2, 0, 0>},
{fft_radix_2<float2, 25, 0, 0, 0>, fft_radix_2<float2, 25, 1, 0, 0>, fft_radix_2<float2, 25, 2, 0, 0>},

};
};


template<> struct TurboFFT_Kernel_Entry<float2, 1, 0>
{
void (*turboFFTArr [26][3])(float2 *, float2 *, float2 *, float2 *, int, int) ={
 {NULL, NULL, NULL},
{fft_radix_2<float2, 1, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 2, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 3, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 4, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 5, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 6, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 7, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 8, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 9, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 10, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 11, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 12, 0, 1, 0>, NULL, NULL},
{fft_radix_2<float2, 13, 0, 1, 0>, fft_radix_2<float2, 13, 1, 1, 0>, NULL},
{fft_radix_2<float2, 14, 0, 1, 0>, fft_radix_2<float2, 14, 1, 1, 0>, NULL},
{fft_radix_2<float2, 15, 0, 1, 0>, fft_radix_2<float2, 15, 1, 1, 0>, NULL},
{fft_radix_2<float2, 16, 0, 1, 0>, fft_radix_2<float2, 16, 1, 1, 0>, NULL},
{fft_radix_2<float2, 17, 0, 1, 0>, fft_radix_2<float2, 17, 1, 1, 0>, NULL},
{fft_radix_2<float2, 18, 0, 1, 0>, fft_radix_2<float2, 18, 1, 1, 0>, NULL},
{fft_radix_2<float2, 19, 0, 1, 0>, fft_radix_2<float2, 19, 1, 1, 0>, NULL},
{fft_radix_2<float2, 20, 0, 1, 0>, fft_radix_2<float2, 20, 1, 1, 0>, NULL},
{fft_radix_2<float2, 21, 0, 1, 0>, fft_radix_2<float2, 21, 1, 1, 0>, NULL},
{fft_radix_2<float2, 22, 0, 1, 0>, fft_radix_2<float2, 22, 1, 1, 0>, NULL},
{fft_radix_2<float2, 23, 0, 1, 0>, fft_radix_2<float2, 23, 1, 1, 0>, fft_radix_2<float2, 23, 2, 1, 0>},
{fft_radix_2<float2, 24, 0, 1, 0>, fft_radix_2<float2, 24, 1, 1, 0>, fft_radix_2<float2, 24, 2, 1, 0>},
{fft_radix_2<float2, 25, 0, 1, 0>, fft_radix_2<float2, 25, 1, 1, 0>, fft_radix_2<float2, 25, 2, 1, 0>},

};
};


template<> struct TurboFFT_Kernel_Entry<float2, 1, 1>
{
void (*turboFFTArr [26][3])(float2 *, float2 *, float2 *, float2 *, int, int) ={
 {NULL, NULL, NULL},
{fft_radix_2<float2, 1, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 2, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 3, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 4, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 5, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 6, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 7, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 8, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 9, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 10, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 11, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 12, 0, 1, 1>, NULL, NULL},
{fft_radix_2<float2, 13, 0, 1, 1>, fft_radix_2<float2, 13, 1, 1, 1>, NULL},
{fft_radix_2<float2, 14, 0, 1, 1>, fft_radix_2<float2, 14, 1, 1, 1>, NULL},
{fft_radix_2<float2, 15, 0, 1, 1>, fft_radix_2<float2, 15, 1, 1, 1>, NULL},
{fft_radix_2<float2, 16, 0, 1, 1>, fft_radix_2<float2, 16, 1, 1, 1>, NULL},
{fft_radix_2<float2, 17, 0, 1, 1>, fft_radix_2<float2, 17, 1, 1, 1>, NULL},
{fft_radix_2<float2, 18, 0, 1, 1>, fft_radix_2<float2, 18, 1, 1, 1>, NULL},
{fft_radix_2<float2, 19, 0, 1, 1>, fft_radix_2<float2, 19, 1, 1, 1>, NULL},
{fft_radix_2<float2, 20, 0, 1, 1>, fft_radix_2<float2, 20, 1, 1, 1>, NULL},
{fft_radix_2<float2, 21, 0, 1, 1>, fft_radix_2<float2, 21, 1, 1, 1>, NULL},
{fft_radix_2<float2, 22, 0, 1, 1>, fft_radix_2<float2, 22, 1, 1, 1>, NULL},
{fft_radix_2<float2, 23, 0, 1, 1>, fft_radix_2<float2, 23, 1, 1, 1>, fft_radix_2<float2, 23, 2, 1, 1>},
{fft_radix_2<float2, 24, 0, 1, 1>, fft_radix_2<float2, 24, 1, 1, 1>, fft_radix_2<float2, 24, 2, 1, 1>},
{fft_radix_2<float2, 25, 0, 1, 1>, fft_radix_2<float2, 25, 1, 1, 1>, fft_radix_2<float2, 25, 2, 1, 1>},

};
};


template<> struct TurboFFT_Kernel_Entry<double2, 0, 0>
{
void (*turboFFTArr [26][3])(double2 *, double2 *, double2 *, double2 *, int, int) ={
 {NULL, NULL, NULL},
{fft_radix_2<double2, 1, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 2, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 3, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 4, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 5, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 6, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 7, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 8, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 9, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 10, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 11, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 12, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 13, 0, 0, 0>, NULL, NULL},
{fft_radix_2<double2, 14, 0, 0, 0>, fft_radix_2<double2, 14, 1, 0, 0>, NULL},
{fft_radix_2<double2, 15, 0, 0, 0>, fft_radix_2<double2, 15, 1, 0, 0>, NULL},
{fft_radix_2<double2, 16, 0, 0, 0>, fft_radix_2<double2, 16, 1, 0, 0>, NULL},
{fft_radix_2<double2, 17, 0, 0, 0>, fft_radix_2<double2, 17, 1, 0, 0>, NULL},
{fft_radix_2<double2, 18, 0, 0, 0>, fft_radix_2<double2, 18, 1, 0, 0>, NULL},
{fft_radix_2<double2, 19, 0, 0, 0>, fft_radix_2<double2, 19, 1, 0, 0>, NULL},
{fft_radix_2<double2, 20, 0, 0, 0>, fft_radix_2<double2, 20, 1, 0, 0>, NULL},
{fft_radix_2<double2, 21, 0, 0, 0>, fft_radix_2<double2, 21, 1, 0, 0>, NULL},
{fft_radix_2<double2, 22, 0, 0, 0>, fft_radix_2<double2, 22, 1, 0, 0>, NULL},
{fft_radix_2<double2, 23, 0, 0, 0>, fft_radix_2<double2, 23, 1, 0, 0>, fft_radix_2<double2, 23, 2, 0, 0>},
{fft_radix_2<double2, 24, 0, 0, 0>, fft_radix_2<double2, 24, 1, 0, 0>, fft_radix_2<double2, 24, 2, 0, 0>},
{fft_radix_2<double2, 25, 0, 0, 0>, fft_radix_2<double2, 25, 1, 0, 0>, fft_radix_2<double2, 25, 2, 0, 0>},

};
};


template<> struct TurboFFT_Kernel_Entry<double2, 1, 0>
{
void (*turboFFTArr [26][3])(double2 *, double2 *, double2 *, double2 *, int, int) ={
 {NULL, NULL, NULL},
{fft_radix_2<double2, 1, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 2, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 3, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 4, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 5, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 6, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 7, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 8, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 9, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 10, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 11, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 12, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 13, 0, 1, 0>, NULL, NULL},
{fft_radix_2<double2, 14, 0, 1, 0>, fft_radix_2<double2, 14, 1, 1, 0>, NULL},
{fft_radix_2<double2, 15, 0, 1, 0>, fft_radix_2<double2, 15, 1, 1, 0>, NULL},
{fft_radix_2<double2, 16, 0, 1, 0>, fft_radix_2<double2, 16, 1, 1, 0>, NULL},
{fft_radix_2<double2, 17, 0, 1, 0>, fft_radix_2<double2, 17, 1, 1, 0>, NULL},
{fft_radix_2<double2, 18, 0, 1, 0>, fft_radix_2<double2, 18, 1, 1, 0>, NULL},
{fft_radix_2<double2, 19, 0, 1, 0>, fft_radix_2<double2, 19, 1, 1, 0>, NULL},
{fft_radix_2<double2, 20, 0, 1, 0>, fft_radix_2<double2, 20, 1, 1, 0>, NULL},
{fft_radix_2<double2, 21, 0, 1, 0>, fft_radix_2<double2, 21, 1, 1, 0>, NULL},
{fft_radix_2<double2, 22, 0, 1, 0>, fft_radix_2<double2, 22, 1, 1, 0>, NULL},
{fft_radix_2<double2, 23, 0, 1, 0>, fft_radix_2<double2, 23, 1, 1, 0>, fft_radix_2<double2, 23, 2, 1, 0>},
{fft_radix_2<double2, 24, 0, 1, 0>, fft_radix_2<double2, 24, 1, 1, 0>, fft_radix_2<double2, 24, 2, 1, 0>},
{fft_radix_2<double2, 25, 0, 1, 0>, fft_radix_2<double2, 25, 1, 1, 0>, fft_radix_2<double2, 25, 2, 1, 0>},

};
};


template<> struct TurboFFT_Kernel_Entry<double2, 1, 1>
{
void (*turboFFTArr [26][3])(double2 *, double2 *, double2 *, double2 *, int, int) ={
 {NULL, NULL, NULL},
{fft_radix_2<double2, 1, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 2, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 3, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 4, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 5, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 6, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 7, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 8, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 9, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 10, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 11, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 12, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 13, 0, 1, 1>, NULL, NULL},
{fft_radix_2<double2, 14, 0, 1, 1>, fft_radix_2<double2, 14, 1, 1, 1>, NULL},
{fft_radix_2<double2, 15, 0, 1, 1>, fft_radix_2<double2, 15, 1, 1, 1>, NULL},
{fft_radix_2<double2, 16, 0, 1, 1>, fft_radix_2<double2, 16, 1, 1, 1>, NULL},
{fft_radix_2<double2, 17, 0, 1, 1>, fft_radix_2<double2, 17, 1, 1, 1>, NULL},
{fft_radix_2<double2, 18, 0, 1, 1>, fft_radix_2<double2, 18, 1, 1, 1>, NULL},
{fft_radix_2<double2, 19, 0, 1, 1>, fft_radix_2<double2, 19, 1, 1, 1>, NULL},
{fft_radix_2<double2, 20, 0, 1, 1>, fft_radix_2<double2, 20, 1, 1, 1>, NULL},
{fft_radix_2<double2, 21, 0, 1, 1>, fft_radix_2<double2, 21, 1, 1, 1>, NULL},
{fft_radix_2<double2, 22, 0, 1, 1>, fft_radix_2<double2, 22, 1, 1, 1>, NULL},
{fft_radix_2<double2, 23, 0, 1, 1>, fft_radix_2<double2, 23, 1, 1, 1>, fft_radix_2<double2, 23, 2, 1, 1>},
{fft_radix_2<double2, 24, 0, 1, 1>, fft_radix_2<double2, 24, 1, 1, 1>, fft_radix_2<double2, 24, 2, 1, 1>},
{fft_radix_2<double2, 25, 0, 1, 1>, fft_radix_2<double2, 25, 1, 1, 1>, fft_radix_2<double2, 25, 2, 1, 1>},

};
};