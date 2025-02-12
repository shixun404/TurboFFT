namespace turbofft{
namespace fft{
namespace threadblock{
    template<
    typename DataType,
    size_t bs,
    typename Tensor,
    bool transform=false
    >
    __global__ fft(DataType *input, DataType *output){
        
        // Shared to Register

        // ToDo: tensor core
        turbofft::fft::thread::fft();


        // Register to Shared

    }
}
}
}