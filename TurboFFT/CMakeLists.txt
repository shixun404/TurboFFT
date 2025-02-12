cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(turbofft LANGUAGES CXX CUDA)

# Set CUDA architecture
# set(CMAKE_CUDA_ARCHITECTURES 80)

# Set allowed values for ARCH_SM (default to 80)
set(ARCH_SM 80 CACHE STRING "Specify GPU architecture (75 or 80).")
set(ALLOWED_ARCH_SM_VALUES "75;80")

# Validate ARCH_SM
list(FIND ALLOWED_ARCH_SM_VALUES "${ARCH_SM}" ARCH_SM_INDEX)
if(ARCH_SM_INDEX EQUAL -1)
    message(FATAL_ERROR "Invalid ARCH_SM value: ${ARCH_SM}. Allowed values are 75 or 80.")
endif()

# Set CUDA architecture based on ARCH_SM
if(ARCH_SM STREQUAL "75")
    set(CMAKE_CUDA_ARCHITECTURES 75)
elseif(ARCH_SM STREQUAL "80")
    set(CMAKE_CUDA_ARCHITECTURES 80)
endif()


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
# target_compile_definitions(turbofft PRIVATE ARCH_SM=80)

# Pass ARCH_SM to the compiler
target_compile_definitions(turbofft PRIVATE ARCH_SM=${ARCH_SM})


find_package(CUDAToolkit)
target_link_libraries(turbofft CUDA::cudart)
target_link_libraries(turbofft CUDA::cufft)
target_link_libraries(turbofft CUDA::cublas)