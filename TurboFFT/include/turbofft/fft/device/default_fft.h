#include "include/turbofft/tensor.h"
namespace turbofft{
namespace fft{
namespace device{
    template<
    
    /// Input Datatype
    typename DataTypeIn,
    
    /// Output Datatype
    typename DataTypeOut,
    
    /// Device-level Tensor Shape (bs, N1, N2, N3, ...) 
    /// Original Problem Size, ToDo: support multi nodes
    typename Tensor DataTensor_Device,

    /// Kernel-level Tensor Shape (bs, N1, N2, N3, ...)
    /// Single Node Problem Size, N1 * N2 * N3 == N_i(device_level), 
    /// ToDo: bs is constrained by GPU memory size
    /// Global Memory: ToDo
    ///     - glb2sh iterator
    ///     - sh2reg iterator & reg2glb iterator
    ///        - Transform for last itr, (bs, N1, N2, N3) --> (bs, N3, N2, N1)
    typename Tensor DataTensor_Kernel,
    
    /// Threadblock-level Tensor Shape (bs, N1, N2, N3, ...)
    /// Threadblock-level plan, N1 * N2 * N3 = N_i(kernel_level)
    /// bs is limited by shared memory size.
    /// Shared Memory: 
    ///     - assign shared memory load, and access
    ///        - sh2reg iterator
    ///        - reg2sh iterator
    ///     - layout ToDo
    ///        - Each transaction, layout is permuted
    ///        - (N1, N2, N3, bs) --> (N2, N3, N1, bs) --> (N3, N2, N1, bs)
    ///        - bs is always the inner dim to avoid bankconflict
    /// ToDo Twiddling Factor
    
    typename Tensor DataTensor_Threadblock,
    /// Warp-level Tensor Shape (bs, N1, N2, N3, ...)
    /// Warp-level plan, 
    /// If using tensor core, 
    ///       N1 * N2 * N3 ... = N_i(threadblock_level) 
    /// Else, 
    ///       No plan is required.
    /// bs = WarpSize, namely 32 by default
    /// 
    
    typename Tensor DataTensor_Warp,
    /// Thread-level Tensor Shape (bs, N1, N2, N3, ...)
    /// Thread-level plan
    /// If using tensor core, 
    ///       No plan is required.
    /// Else, 
    ///       N1 * N2 * N3 ... = N_i(threadblock_level) 
    /// bs = N_i(thredblock_level) / (N1 * N2 * N3 * ...)
    typename Tensor DataTensor_Thread,

    /// Global to register iterator
    typename glb2reg_itr,

    /// Register to global iterator
    typename reg2glb_itr,

    /// Register to shared iterator, global
    typename reg2sh_glb_itr,
    
    /// Shared to register iterator, global
    typename reg2sh_glb_itr,

    /// Register to shared iterator
    typename reg2sh_itr,
    
    /// Shared to register iterator
    typename reg2sh_itr,

    /// Operator
    typename Operator
    > struct DefaultFFT{
        using FFTKernel = 
    }
}
}
}