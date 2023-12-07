#include "include/turbofft/tensor.h"
#include <stdio.h>
int main(){
    turbofft::fft::Tensor<float, 2, 2, 2, 2> a;
    printf("Size=%d\n", a.CalculateTotalSize());
    return 0;
}