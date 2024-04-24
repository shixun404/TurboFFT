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
#include "code_gen/generated/double2/fft_radix_2_logN_21_upload_2.cuh"
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
#include "code_gen/generated/double2/fft_radix_2_logN_26_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_26_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_26_upload_2.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_27_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_27_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_27_upload_2.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_28_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_28_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_28_upload_2.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_29_upload_0.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_29_upload_1.cuh"
#include "code_gen/generated/double2/fft_radix_2_logN_29_upload_2.cuh"


template <>
struct TurboFFT_Kernel_Entry<double2> {
void (*turboFFTArr[30][3])(double2 *, double2 *, double2 *, double2*, int, int) = {
    {NULL, NULL, NULL},
    {fft_radix_2_logN<double2, 1, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 2, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 3, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 4, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 5, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 6, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 7, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 8, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 9, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 10, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 11, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 12, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 13, 0>, NULL, NULL},
    {fft_radix_2_logN<double2, 14, 0>, fft_radix_2_logN<double2, 14, 1>, NULL},
    {fft_radix_2_logN<double2, 15, 0>, fft_radix_2_logN<double2, 15, 1>, NULL},
    {fft_radix_2_logN<double2, 16, 0>, fft_radix_2_logN<double2, 16, 1>, NULL},
    {fft_radix_2_logN<double2, 17, 0>, fft_radix_2_logN<double2, 17, 1>, NULL},
    {fft_radix_2_logN<double2, 18, 0>, fft_radix_2_logN<double2, 18, 1>, NULL},
    {fft_radix_2_logN<double2, 19, 0>, fft_radix_2_logN<double2, 19, 1>, NULL},
    {fft_radix_2_logN<double2, 20, 0>, fft_radix_2_logN<double2, 20, 1>, NULL},
    {fft_radix_2_logN<double2, 21, 0>, fft_radix_2_logN<double2, 21, 1>, NULL},
    {fft_radix_2_logN<double2, 22, 0>, fft_radix_2_logN<double2, 22, 1>, fft_radix_2_logN<double2, 22, 2>},
    {fft_radix_2_logN<double2, 23, 0>, fft_radix_2_logN<double2, 23, 1>, fft_radix_2_logN<double2, 23, 2>},
    {fft_radix_2_logN<double2, 24, 0>, fft_radix_2_logN<double2, 24, 1>, fft_radix_2_logN<double2, 24, 2>},
    {fft_radix_2_logN<double2, 25, 0>, fft_radix_2_logN<double2, 25, 1>, fft_radix_2_logN<double2, 25, 2>}
};
}