#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <cuda_fp16.h>
#include <cuda_bf16.h>


// #define CUDA_CALLER(call) do{\
//   cudaError_t cuda_ret = (call);\
//   if(cuda_ret != cudaSuccess){\
//     printf("CUDA Error at line %d in file %s\n", __LINE__, __FILE__);\
//     printf("  Error message: %s\n", cudaGetErrorString(cuda_ret));\
//     printf("  In the function call %s\n", #call);\
//     exit(1);\
//   }\
// }while(0)

// // #define checkCudaErrors(err) __checkCudaErrors(err, __FILE__, __LINE__)

// // // These are the inline versions for all of the SDK helper functions
// // inline void __checkCudaErrors(CUresult err, const char *file, const int line) {
// //   if (CUDA_SUCCESS != err) {
// //     const char *errorStr = NULL;
// //     cuGetErrorString(err, &errorStr);
// //     fprintf(stderr,
// //             "checkCudaErrors() Driver API error = %04d \"%s\" from file <%s>, "
// //             "line %i.\n",
// //             err, errorStr, file, line);
// //     exit(EXIT_FAILURE);
// //   }
// // }