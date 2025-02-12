#include "include/turbofft/tensor.h"
#include "include/turbofft/fft/threadblock/tensor.h"
namespace turbofft{
namespace fft{
namespace kernel{
    template<
    typename DataType,
    size_t bs,
    typename Tensor,
    bool transform=false
    >
    __global__ fft(DataType *input, DataType *output){
        
        
        Tensor a(input, );
        
        /// Global to Shared


        /// Threadblock FFT
        turbofft::fft::threadblock::fft();

        /// Shared to Global

        

    }
}
}
}