namespace turbofft {

namespace fft {

namespace kernel {

////////////////////////////////////////////////////////////////////////////////

template<
    /// Input Datatype: Fload, Double, Complex, ComplexDouble
    typename DataTypeIn,
    /// Output Datatype: Fload, Double, Complex, ComplexDouble
    typename DataTypeOut,
    /// N-point FFT in Kernel level
    int n_kernel,
    /// n-point FFT in threadblock level
    int n_threadblock,
    /// n-point FFT in warp level
    int n_warp,
    /// n-point FFT in thread level
    int n_thread,
    /// Batch size for kernel-level
    int bs_kernel,
    /// Batch Size for threadblock-level
    int bs_threadbock,
    /// Batch Size for warp-level
    int bs_warp,
    /// Batch Size for thread-level
    int bs_thread,
    /// Stride n direction in global memory
    int stride_global_n,
    /// Stride bs direction in global memory
    int stride_global_bs,
    /// Stride n direction in shared memory
    int stride_shared_n,
    /// Stride bs direction in shared memory
    int stride_shared_bs,
    /// Stride between continuous signals in a batch
    int stride_batch,
    
    /// Shared memory padding
    int shmem_padding,
    /// Operator class tag
    typename OperatorClass,
    /// Tag indicating architecture to tune for
    typename ArchTag,
    /// Transform Output
    bool Transform,
>
struct DefaultFFT;

////////////////////////////////////////////////////////////////////////////////



}

}

}