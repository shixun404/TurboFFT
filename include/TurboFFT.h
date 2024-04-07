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
#include "utils/abft.h"
#include "utils/CommandLineParser.h"

template<typename DataType, int N, int dim_id> void __global__ fft_radix_2(DataType *, DataType *, DataType *, DataType*, int, int){

};

template <typename DataType>
struct TurboFFT_Kernel_Entry {
void (*turboFFTArr[30][3])(DataType *, DataType *, DataType *, DataType*, int, int);
};

template <typename DataType>
void test_turbofft( DataType* input_d, DataType* output_d, DataType* output_turbofft,
                    DataType* twiddle_d, DataType* checksum, std::vector<long long int> param, 
                    long long int bs, int ntest, auto alpha);