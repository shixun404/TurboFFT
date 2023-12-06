#include "include/turbofft/tensor.h"

int main(){
    int dim = 2;
    int dims[2] = {2, 2};
    turbofft::fft::Tensor<float, dim, dims>;
    return 0;
}