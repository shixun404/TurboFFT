template<typename DataType, int N, int dim_id> void __global__ fft_radix_2(DataType *, DataType *, DataType *, DataType*, int, int);
extern __shared__ float ext_shared[];