cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(turbofft LANGUAGES CXX CUDA)

# Set CUDA architecture
set(CMAKE_CUDA_ARCHITECTURES 80)

# Add -O3 optimization flag for CUDA
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -lineinfo")

# Find CUDA package
find_package(CUDA REQUIRED)

# Include the TurboFFT directory
include_directories("include/")

include_directories("~/cuda-samples/Common/")

# Specify the executable and its source files
add_executable(turbofft main.cu)

# Link the CUDA libraries
# set(CUDA_TOOLKIT_ROOT_DIR "/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/12.0")
# # set(CUDA_INCLUDE_DIRS "/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/12.0/include")
# set(CUDA_CUDART_LIBRARY "/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/12.0/lib64")
find_package(CUDAToolkit)
target_link_libraries(turbofft CUDA::cudart)
target_link_libraries(turbofft CUDA::cufft)
target_link_libraries(turbofft CUDA::cublas)