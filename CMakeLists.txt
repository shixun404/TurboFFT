cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(cuda_project LANGUAGES CXX CUDA)

# Set CUDA architecture
set(CMAKE_CUDA_ARCHITECTURES 80)

# Find CUDA package
find_package(CUDA REQUIRED)

# Include the TurboFFT directory
include_directories("include")

# Specify the executable and its source files
add_executable(cuda_project test.cu)

# Link the CUDA libraries
target_link_libraries(cuda_project CUDA::cudart)