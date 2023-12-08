#include <stdio.h>
#include <cuda_runtime.h> 
#include <cufftXt.h>
#include <cuda_fp16.h>
#include <cuda_bf16.h>

#include "utils/utils.h"
#include "turbofft/fft/thread/fft.h"
#include "profiler/cufft/cufft.h"
#include "utils/compareData.h"
#include "utils/printData.h"
#include "utils/initializeData.h"