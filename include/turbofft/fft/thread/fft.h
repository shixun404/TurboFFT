#include "turbofft/tensor.h"
#include "stdio.h"
namespace turbofft{
namespace fft{
namespace thread{
    template<
    typename DataType,
    typename Tensor
    >
    // template<typename DataType>
    __global__ void fft(DataType *input, DataType *output){
        DataType tmp[2];

        
        tmp[0].x = input[0].x + input[1].x;
        tmp[0].y = input[0].y + input[1].y;
              
        tmp[1].x = input[0].x - input[1].x;
        tmp[1].y = input[0].y - input[1].y;

        output[0] = tmp[0];
        output[1] = tmp[1];

    }
    
////////////////////////////////////////////////////////////////////////////////
/// Partial Specialization for 2-point    
    // template<
    // typename DataType
    // >
    // __global__ void fft<DataType, turbofft::Tensor<DataType, 1, 2>>(DataType *input, DataType *output){
        
    //     DataType tmp[2];
        
    //     tmp[0].x = input[0].x + input[1].x;
    //     tmp[0].y = input[0].y + input[1].y;
        
    //     tmp[1].x = input[0].x - input[1].x;
    //     tmp[1].y = input[0].y - input[1].y;

    //     output[0] = tmp[0];
    //     output[1] = tmp[1];
    // }
// ////////////////////////////////////////////////////////////////////////////////
// /// Partial Specialization for 4-point
//     template<
//     typename DataType
//     >
//     __global__ fft<DataType, bs, Tensor<DataType, bs, 2, 2>>(DataType *input, DataType *output){

//     }
// ////////////////////////////////////////////////////////////////////////////////
// /// Partial Specialization for 8-point
//     template<
//     typename DataType
//     >
//     __global__ fft<DataType, bs, Tensor<DataType, bs, 2, 2, 2>>
//     (DataType *input, DataType *output){

//     }
// ////////////////////////////////////////////////////////////////////////////////
// /// Partial Specialization for 16-point
//     template<
//     typename DataType
//     >
//     __global__ fft<DataType, bs, Tensor<DataType, bs, 2, 2, 2, 2>>
//     (DataType *input, DataType *output){

//     }
// ////////////////////////////////////////////////////////////////////////////////
// /// Partial Specialization for 32-point
//     template<
//     typename DataType
//     >
//     __global__ fft<DataType, bs, Tensor<DataType, bs, 2, 2, 2, 2, 2>>
//     (DataType *input, DataType *output){

//     }
// ////////////////////////////////////////////////////////////////////////////////
// /// Partial Specialization for 64-point
//     template<
//     typename DataType
//     >
//     __global__ fft<DataType, bs, Tensor<DataType, bs, 2, 2, 2, 2, 2, 2>>
//     (DataType *input, DataType *output){

//     }
}
}
}