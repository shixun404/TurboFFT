#define CUDA_CALLER(call) do{\
  cudaError_t cuda_ret = (call);\
  if(cuda_ret != cudaSuccess){\
    printf("CUDA Error at line %d in file %s\n", __LINE__, __FILE__);\
    printf("  Error message: %s\n", cudaGetErrorString(cuda_ret));\
    printf("  In the function call %s\n", #call);\
    exit(1);\
  }\
}while(0)