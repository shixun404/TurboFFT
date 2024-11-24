cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(turbofft LANGUAGES CXX CUDA)

# Set CUDA architecture
set(CMAKE_CUDA_ARCHITECTURES 80)

# Add -O3 optimization flag for CUDA
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -lineinfo")

# Find CUDA package
find_package(CUDA REQUIRED)

# Include the TurboFFT directory
include_directories("include/")

# Include CUDA Samples Common via environment variable
if(DEFINED ENV{CUDA_SAMPLES_COMMON})
    include_directories("$ENV{CUDA_SAMPLES_COMMON}")
else()
    message(FATAL_ERROR "CUDA Samples Common directory not set. Define CUDA_SAMPLES_COMMON environment variable.")
endif()

# Specify the executable and its source files
add_executable(turbofft main.cu)
target_compile_definitions(turbofft PRIVATE ARCH_SM=80)

find_package(CUDAToolkit)
target_link_libraries(turbofft CUDA::cudart)
target_link_libraries(turbofft CUDA::cufft)
target_link_libraries(turbofft CUDA::cublas)