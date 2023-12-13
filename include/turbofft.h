#include <stdio.h>
#include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cuda_fp16.h>
#include <cuda_bf16.h>

#include "utils/utils.h"
#include "turbofft/fft/thread/fft.h"
#include "cufft/cufft.h"
#include "cufft/cufft_ft.h"
#include "utils/compareData.h"
#include "utils/printData.h"
#include "utils/initializeData.h"