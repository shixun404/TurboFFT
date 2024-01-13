#include <stdio.h>
#include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cuda_fp16.h>
#include <cuda_bf16.h>
#include <fstream>

#include "utils/utils.h"
#include "turbofft/fft/thread/fft.h"
#include "cufft/cufft.h"
#include "cufft/cufft_ft.h"
#include "utils/compareData.h"
#include "utils/printData.h"
#include "utils/initializeData.h"
#include "utils/readCSV.h"

#include "code_gen/generated/fft_radix_2_logN_1_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_2_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_3_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_4_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_5_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_6_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_7_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_8_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_9_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_10_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_11_upload_0.cuh"
#include "code_gen/generated/fft_radix_2_logN_12_upload_0.cuh"

void (*turboFFTArr[13][3])(double2 *, double2 *, double2 *, int) = {
    {NULL, NULL, NULL},
    {fft_radix_2_logN_1_dim_0, NULL, NULL},
    {fft_radix_2_logN_2_dim_0, NULL, NULL},
    {fft_radix_2_logN_3_dim_0, NULL, NULL},
    {fft_radix_2_logN_4_dim_0, NULL, NULL},
    {fft_radix_2_logN_5_dim_0, NULL, NULL},
    {fft_radix_2_logN_6_dim_0, NULL, NULL},
    {fft_radix_2_logN_7_dim_0, NULL, NULL},
    {fft_radix_2_logN_8_dim_0, NULL, NULL},
    {fft_radix_2_logN_9_dim_0, NULL, NULL},
    {fft_radix_2_logN_10_dim_0, NULL, NULL},
    {fft_radix_2_logN_11_dim_0, NULL, NULL},
    {fft_radix_2_logN_12_dim_0, NULL, NULL}
};
