#include "turbofft/tensor.h"
#include "turbofft/constants.h"
#include "turbofft/macro_ops.h"
#include "stdio.h"
namespace turbofft{
namespace fft{
namespace thread{
    extern __shared__ double shared[];

    template<
    typename DataType,
    int logN
    >
    __global__ void fft(DataType *input, DataType *output);

    // template<typename DataType, 1>
    // __global__ void fft<DataType, 1>(DataType *input, DataType *output){
    //     DataType tmp[2];

        
    //     tmp[0].x = input[0].x + input[1].x;
    //     tmp[0].y = input[0].y + input[1].y;
              
    //     tmp[1].x = input[0].x - input[1].x;
    //     tmp[1].y = input[0].y - input[1].y;

    //     output[0] = tmp[0];
    //     output[1] = tmp[1];

    // }

    
    
    template<
    typename DataType,
    int logN
    >
    __global__ void  fft(DataType* inputs, DataType* outputs) {
    
        int tid = threadIdx.x + threadIdx.y * blockDim.x;
        DataType temp_0;
        DataType temp_1;
        DataType temp_2;
        DataType temp_3;
        DataType temp_4;
        DataType temp_5;
        DataType temp_6;
        DataType temp_7;
        
        DataType* sdata = (DataType*)shared;
        int tx = threadIdx.x;
        int ty = threadIdx.y;
        int bx = blockIdx.x;
        int N = 256;
        int __id[8];
        DataType tmp;
        DataType tmp_angle, tmp_angle_rot;
        int j;
        int k;
        int tmp_id;
        int n = 1, n_global = 1;
        DataType tmp_angle_bk;
        temp_0 = inputs[(ty + 0 * 32) + (tx + bx * 1) * 256];
        temp_1 = inputs[(ty + 1 * 32) + (tx + bx * 1) * 256];
        temp_2 = inputs[(ty + 2 * 32) + (tx + bx * 1) * 256];
        temp_3 = inputs[(ty + 3 * 32) + (tx + bx * 1) * 256];
        temp_4 = inputs[(ty + 4 * 32) + (tx + bx * 1) * 256];
        temp_5 = inputs[(ty + 5 * 32) + (tx + bx * 1) * 256];
        temp_6 = inputs[(ty + 6 * 32) + (tx + bx * 1) * 256];
        temp_7 = inputs[(ty + 7 * 32) + (tx + bx * 1) * 256];
        
        __id[0] = 0 + ty;
        __id[1] = 32 + ty;
        __id[2] = 64 + ty;
        __id[3] = 96 + ty;
        __id[4] = 128 + ty;
        __id[5] = 160 + ty;
        __id[6] = 192 + ty;
        __id[7] = 224 + ty;
        
        j = 1;
        k = 4 % 1;
        MY_ANGLE2COMPLEX((double)(j * k) * -3.141592653589793f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
            tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_4, tmp_angle, tmp);
            temp_4 = tmp;
            
            MY_MUL(temp_5, tmp_angle, tmp);
            temp_5 = tmp;
            
            MY_MUL(temp_6, tmp_angle, tmp);
            temp_6 = tmp;
            
            MY_MUL(temp_7, tmp_angle, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_4, temp_0);
            MY_SUB(tmp, temp_4, temp_4);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[4] = tmp_id + 1;
            
            tmp = temp_1;
            MY_ADD(tmp, temp_5, temp_1);
            MY_SUB(tmp, temp_5, temp_5);
            
            tmp_id = __id[1];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[1] = tmp_id;
            __id[5] = tmp_id + 1;
            
            tmp = temp_2;
            MY_ADD(tmp, temp_6, temp_2);
            MY_SUB(tmp, temp_6, temp_6);
            
            tmp_id = __id[2];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[2] = tmp_id;
            __id[6] = tmp_id + 1;
            
            tmp = temp_3;
            MY_ADD(tmp, temp_7, temp_3);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[3];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[3] = tmp_id;
            __id[7] = tmp_id + 1;
            
        n_global *= 2;
        
        j = 1;
        k = 4 % 2;
        MY_ANGLE2COMPLEX((double)(j * k) * -1.5707963267948966f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
            tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_2, tmp_angle, tmp);
            temp_2 = tmp;
            
            MY_MUL(temp_6, tmp_angle_rot, tmp);
            temp_6 = tmp;
            
            MY_MUL(temp_3, tmp_angle, tmp);
            temp_3 = tmp;
            
            MY_MUL(temp_7, tmp_angle_rot, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_2, temp_0);
            MY_SUB(tmp, temp_2, temp_2);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[2] = tmp_id + 2;
            
            tmp = temp_4;
            MY_ADD(tmp, temp_6, temp_4);
            MY_SUB(tmp, temp_6, temp_6);
            
            tmp_id = __id[4];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[4] = tmp_id;
            __id[6] = tmp_id + 2;
            
            tmp = temp_1;
            MY_ADD(tmp, temp_3, temp_1);
            MY_SUB(tmp, temp_3, temp_3);
            
            tmp_id = __id[1];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[1] = tmp_id;
            __id[3] = tmp_id + 2;
            
            tmp = temp_5;
            MY_ADD(tmp, temp_7, temp_5);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[5];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[5] = tmp_id;
            __id[7] = tmp_id + 2;
            
        n_global *= 2;
        
        j = 1;
        k = 4 % 4;
        MY_ANGLE2COMPLEX((double)(j * k) * -0.7853981633974483f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_1, tmp_angle, tmp);
            temp_1 = tmp;
            
            MY_MUL(temp_3, tmp_angle_rot, tmp);
            temp_3 = tmp;
            
            tmp_angle_rot.x = 0.7071067811865476f;
            tmp_angle_rot.y = -0.7071067811865475f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_5, tmp_angle, tmp);
            temp_5 = tmp;
            
            MY_MUL(temp_7, tmp_angle_rot, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_1, temp_0);
            MY_SUB(tmp, temp_1, temp_1);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[1] = tmp_id + 4;
            
            tmp = temp_4;
            MY_ADD(tmp, temp_5, temp_4);
            MY_SUB(tmp, temp_5, temp_5);
            
            tmp_id = __id[4];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[4] = tmp_id;
            __id[5] = tmp_id + 4;
            
            tmp = temp_2;
            MY_ADD(tmp, temp_3, temp_2);
            MY_SUB(tmp, temp_3, temp_3);
            
            tmp_id = __id[2];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[2] = tmp_id;
            __id[3] = tmp_id + 4;
            
            tmp = temp_6;
            MY_ADD(tmp, temp_7, temp_6);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[6];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[6] = tmp_id;
            __id[7] = tmp_id + 4;
            
        n_global *= 2;
        
        __syncthreads();
        
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 0) / (double)(256), tmp_angle);
    MY_MUL(temp_0, tmp_angle, tmp);
    temp_0 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 1) / (double)(256), tmp_angle);
    MY_MUL(temp_4, tmp_angle, tmp);
    temp_4 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 2) / (double)(256), tmp_angle);
    MY_MUL(temp_2, tmp_angle, tmp);
    temp_2 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 3) / (double)(256), tmp_angle);
    MY_MUL(temp_6, tmp_angle, tmp);
    temp_6 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 4) / (double)(256), tmp_angle);
    MY_MUL(temp_1, tmp_angle, tmp);
    temp_1 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 5) / (double)(256), tmp_angle);
    MY_MUL(temp_5, tmp_angle, tmp);
    temp_5 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 6) / (double)(256), tmp_angle);
    MY_MUL(temp_3, tmp_angle, tmp);
    temp_3 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 1) * 7) / (double)(256), tmp_angle);
    MY_MUL(temp_7, tmp_angle, tmp);
    temp_7 = tmp;
    
        sdata[((tx + 1 * __id[0]) / 16) * 17 + 
        ((tx + 1 * __id[0]) % 16)] = temp_0;
        
        sdata[((tx + 1 * __id[4]) / 16) * 17 + 
        ((tx + 1 * __id[4]) % 16)] = temp_4;
        
        sdata[((tx + 1 * __id[2]) / 16) * 17 + 
        ((tx + 1 * __id[2]) % 16)] = temp_2;
        
        sdata[((tx + 1 * __id[6]) / 16) * 17 + 
        ((tx + 1 * __id[6]) % 16)] = temp_6;
        
        sdata[((tx + 1 * __id[1]) / 16) * 17 + 
        ((tx + 1 * __id[1]) % 16)] = temp_1;
        
        sdata[((tx + 1 * __id[5]) / 16) * 17 + 
        ((tx + 1 * __id[5]) % 16)] = temp_5;
        
        sdata[((tx + 1 * __id[3]) / 16) * 17 + 
        ((tx + 1 * __id[3]) % 16)] = temp_3;
        
        sdata[((tx + 1 * __id[7]) / 16) * 17 + 
        ((tx + 1 * __id[7]) % 16)] = temp_7;
        
        __syncthreads();		
        
        temp_0 = sdata[((tx + 1 * (0 + ty)) / 16) * 17 +
                            ((tx + 1 * (0 + ty)) % 16)];
        __id[0] = ty + 0;
        
        temp_1 = sdata[((tx + 1 * (32 + ty)) / 16) * 17 +
                            ((tx + 1 * (32 + ty)) % 16)];
        __id[1] = ty + 32;
        
        temp_2 = sdata[((tx + 1 * (64 + ty)) / 16) * 17 +
                            ((tx + 1 * (64 + ty)) % 16)];
        __id[2] = ty + 64;
        
        temp_3 = sdata[((tx + 1 * (96 + ty)) / 16) * 17 +
                            ((tx + 1 * (96 + ty)) % 16)];
        __id[3] = ty + 96;
        
        temp_4 = sdata[((tx + 1 * (128 + ty)) / 16) * 17 +
                            ((tx + 1 * (128 + ty)) % 16)];
        __id[4] = ty + 128;
        
        temp_5 = sdata[((tx + 1 * (160 + ty)) / 16) * 17 +
                            ((tx + 1 * (160 + ty)) % 16)];
        __id[5] = ty + 160;
        
        temp_6 = sdata[((tx + 1 * (192 + ty)) / 16) * 17 +
                            ((tx + 1 * (192 + ty)) % 16)];
        __id[6] = ty + 192;
        
        temp_7 = sdata[((tx + 1 * (224 + ty)) / 16) * 17 +
                            ((tx + 1 * (224 + ty)) % 16)];
        __id[7] = ty + 224;
        
        j = 1;
        k = 4 % 1;
        MY_ANGLE2COMPLEX((double)(j * k) * -3.141592653589793f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_4, tmp_angle, tmp);
            temp_4 = tmp;
            
            MY_MUL(temp_5, tmp_angle, tmp);
            temp_5 = tmp;
            
            MY_MUL(temp_6, tmp_angle, tmp);
            temp_6 = tmp;
            
            MY_MUL(temp_7, tmp_angle, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_4, temp_0);
            MY_SUB(tmp, temp_4, temp_4);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[4] = tmp_id + 8;
            
            tmp = temp_1;
            MY_ADD(tmp, temp_5, temp_1);
            MY_SUB(tmp, temp_5, temp_5);
            
            tmp_id = __id[1];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[1] = tmp_id;
            __id[5] = tmp_id + 8;
            
            tmp = temp_2;
            MY_ADD(tmp, temp_6, temp_2);
            MY_SUB(tmp, temp_6, temp_6);
            
            tmp_id = __id[2];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[2] = tmp_id;
            __id[6] = tmp_id + 8;
            
            tmp = temp_3;
            MY_ADD(tmp, temp_7, temp_3);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[3];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[3] = tmp_id;
            __id[7] = tmp_id + 8;
            
        n_global *= 2;
        
        j = 1;
        k = 4 % 2;
        MY_ANGLE2COMPLEX((double)(j * k) * -1.5707963267948966f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_2, tmp_angle, tmp);
            temp_2 = tmp;
            
            MY_MUL(temp_6, tmp_angle_rot, tmp);
            temp_6 = tmp;
            
            MY_MUL(temp_3, tmp_angle, tmp);
            temp_3 = tmp;
            
            MY_MUL(temp_7, tmp_angle_rot, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_2, temp_0);
            MY_SUB(tmp, temp_2, temp_2);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[2] = tmp_id + 16;
            
            tmp = temp_4;
            MY_ADD(tmp, temp_6, temp_4);
            MY_SUB(tmp, temp_6, temp_6);
            
            tmp_id = __id[4];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[4] = tmp_id;
            __id[6] = tmp_id + 16;
            
            tmp = temp_1;
            MY_ADD(tmp, temp_3, temp_1);
            MY_SUB(tmp, temp_3, temp_3);
            
            tmp_id = __id[1];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[1] = tmp_id;
            __id[3] = tmp_id + 16;
            
            tmp = temp_5;
            MY_ADD(tmp, temp_7, temp_5);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[5];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[5] = tmp_id;
            __id[7] = tmp_id + 16;
            
        n_global *= 2;
        
        j = 1;
        k = 4 % 4;
        MY_ANGLE2COMPLEX((double)(j * k) * -0.7853981633974483f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_1, tmp_angle, tmp);
            temp_1 = tmp;
            
            MY_MUL(temp_3, tmp_angle_rot, tmp);
            temp_3 = tmp;
            
            tmp_angle_rot.x = 0.7071067811865476f;
            tmp_angle_rot.y = -0.7071067811865475f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_5, tmp_angle, tmp);
            temp_5 = tmp;
            
            MY_MUL(temp_7, tmp_angle_rot, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_1, temp_0);
            MY_SUB(tmp, temp_1, temp_1);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[1] = tmp_id + 32;
            
            tmp = temp_4;
            MY_ADD(tmp, temp_5, temp_4);
            MY_SUB(tmp, temp_5, temp_5);
            
            tmp_id = __id[4];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[4] = tmp_id;
            __id[5] = tmp_id + 32;
            
            tmp = temp_2;
            MY_ADD(tmp, temp_3, temp_2);
            MY_SUB(tmp, temp_3, temp_3);
            
            tmp_id = __id[2];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[2] = tmp_id;
            __id[3] = tmp_id + 32;
            
            tmp = temp_6;
            MY_ADD(tmp, temp_7, temp_6);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[6];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[6] = tmp_id;
            __id[7] = tmp_id + 32;
            
        n_global *= 2;
        
        __syncthreads();
        
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 0) / (double)(32.0), tmp_angle);
    MY_MUL(temp_0, tmp_angle, tmp);
    temp_0 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 1) / (double)(32.0), tmp_angle);
    MY_MUL(temp_4, tmp_angle, tmp);
    temp_4 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 2) / (double)(32.0), tmp_angle);
    MY_MUL(temp_2, tmp_angle, tmp);
    temp_2 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 3) / (double)(32.0), tmp_angle);
    MY_MUL(temp_6, tmp_angle, tmp);
    temp_6 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 4) / (double)(32.0), tmp_angle);
    MY_MUL(temp_1, tmp_angle, tmp);
    temp_1 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 5) / (double)(32.0), tmp_angle);
    MY_MUL(temp_5, tmp_angle, tmp);
    temp_5 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 6) / (double)(32.0), tmp_angle);
    MY_MUL(temp_3, tmp_angle, tmp);
    temp_3 = tmp;
    
    MY_ANGLE2COMPLEX((double)(-M_PI * 2 * ((ty) / 8) * 7) / (double)(32.0), tmp_angle);
    MY_MUL(temp_7, tmp_angle, tmp);
    temp_7 = tmp;
    
        sdata[((tx + 1 * __id[0]) / 16) * 17 + 
        ((tx + 1 * __id[0]) % 16)] = temp_0;
        
        sdata[((tx + 1 * __id[4]) / 16) * 17 + 
        ((tx + 1 * __id[4]) % 16)] = temp_4;
        
        sdata[((tx + 1 * __id[2]) / 16) * 17 + 
        ((tx + 1 * __id[2]) % 16)] = temp_2;
        
        sdata[((tx + 1 * __id[6]) / 16) * 17 + 
        ((tx + 1 * __id[6]) % 16)] = temp_6;
        
        sdata[((tx + 1 * __id[1]) / 16) * 17 + 
        ((tx + 1 * __id[1]) % 16)] = temp_1;
        
        sdata[((tx + 1 * __id[5]) / 16) * 17 + 
        ((tx + 1 * __id[5]) % 16)] = temp_5;
        
        sdata[((tx + 1 * __id[3]) / 16) * 17 + 
        ((tx + 1 * __id[3]) % 16)] = temp_3;
        
        sdata[((tx + 1 * __id[7]) / 16) * 17 + 
        ((tx + 1 * __id[7]) % 16)] = temp_7;
        
        __syncthreads();		
        
        temp_0 = sdata[((tx + 1 * (0 + ty)) / 16) * 17 +
                            ((tx + 1 * (0 + ty)) % 16)];
        __id[0] = ty + 0;
        
        temp_1 = sdata[((tx + 1 * (32 + ty)) / 16) * 17 +
                            ((tx + 1 * (32 + ty)) % 16)];
        __id[1] = ty + 32;
        
        temp_2 = sdata[((tx + 1 * (64 + ty)) / 16) * 17 +
                            ((tx + 1 * (64 + ty)) % 16)];
        __id[2] = ty + 64;
        
        temp_3 = sdata[((tx + 1 * (96 + ty)) / 16) * 17 +
                            ((tx + 1 * (96 + ty)) % 16)];
        __id[3] = ty + 96;
        
        temp_4 = sdata[((tx + 1 * (128 + ty)) / 16) * 17 +
                            ((tx + 1 * (128 + ty)) % 16)];
        __id[4] = ty + 128;
        
        temp_5 = sdata[((tx + 1 * (160 + ty)) / 16) * 17 +
                            ((tx + 1 * (160 + ty)) % 16)];
        __id[5] = ty + 160;
        
        temp_6 = sdata[((tx + 1 * (192 + ty)) / 16) * 17 +
                            ((tx + 1 * (192 + ty)) % 16)];
        __id[6] = ty + 192;
        
        temp_7 = sdata[((tx + 1 * (224 + ty)) / 16) * 17 +
                            ((tx + 1 * (224 + ty)) % 16)];
        __id[7] = ty + 224;
        
        j = 1;
        k = 2 % 1;
        MY_ANGLE2COMPLEX((double)(j * k) * -3.141592653589793f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_4, tmp_angle, tmp);
            temp_4 = tmp;
            
            MY_MUL(temp_6, tmp_angle, tmp);
            temp_6 = tmp;
            
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_5, tmp_angle, tmp);
            temp_5 = tmp;
            
            MY_MUL(temp_7, tmp_angle, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_4, temp_0);
            MY_SUB(tmp, temp_4, temp_4);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[4] = tmp_id + 64;
            
            tmp = temp_2;
            MY_ADD(tmp, temp_6, temp_2);
            MY_SUB(tmp, temp_6, temp_6);
            
            tmp_id = __id[2];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[2] = tmp_id;
            __id[6] = tmp_id + 64;
            
            tmp = temp_1;
            MY_ADD(tmp, temp_5, temp_1);
            MY_SUB(tmp, temp_5, temp_5);
            
            tmp_id = __id[1];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[1] = tmp_id;
            __id[5] = tmp_id + 64;
            
            tmp = temp_3;
            MY_ADD(tmp, temp_7, temp_3);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[3];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[3] = tmp_id;
            __id[7] = tmp_id + 64;
            
        n_global *= 2;
        
        j = 1;
        k = 2 % 2;
        MY_ANGLE2COMPLEX((double)(j * k) * -1.5707963267948966f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_2, tmp_angle, tmp);
            temp_2 = tmp;
            
            MY_MUL(temp_6, tmp_angle_rot, tmp);
            temp_6 = tmp;
            
                    tmp_angle = tmp_angle_bk;
    
            tmp_angle_rot.x = 1.0f;
            tmp_angle_rot.y = 0.0f;
            MY_MUL(tmp_angle, tmp_angle_rot, tmp);
            tmp_angle = tmp;
            tmp_angle_rot.x = tmp_angle.y;
            tmp_angle_rot.y = -tmp_angle.x;
            
            MY_MUL(temp_3, tmp_angle, tmp);
            temp_3 = tmp;
            
            MY_MUL(temp_7, tmp_angle_rot, tmp);
            temp_7 = tmp;
            
            tmp = temp_0;
            MY_ADD(tmp, temp_2, temp_0);
            MY_SUB(tmp, temp_2, temp_2);
            
            tmp_id = __id[0];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[0] = tmp_id;
            __id[2] = tmp_id + 128;
            
            tmp = temp_4;
            MY_ADD(tmp, temp_6, temp_4);
            MY_SUB(tmp, temp_6, temp_6);
            
            tmp_id = __id[4];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[4] = tmp_id;
            __id[6] = tmp_id + 128;
            
            tmp = temp_1;
            MY_ADD(tmp, temp_3, temp_1);
            MY_SUB(tmp, temp_3, temp_3);
            
            tmp_id = __id[1];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[1] = tmp_id;
            __id[3] = tmp_id + 128;
            
            tmp = temp_5;
            MY_ADD(tmp, temp_7, temp_5);
            MY_SUB(tmp, temp_7, temp_7);
            
            tmp_id = __id[5];
            tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
            __id[5] = tmp_id;
            __id[7] = tmp_id + 128;
            
        n_global *= 2;
        
        n_global *= 2;
        
        outputs[(tx + bx * 1) * 256 +  __id[0]] = temp_0;
        
        outputs[(tx + bx * 1) * 256 +  __id[1]] = temp_1;
        
        outputs[(tx + bx * 1) * 256 +  __id[4]] = temp_4;
        
        outputs[(tx + bx * 1) * 256 +  __id[5]] = temp_5;
        
        outputs[(tx + bx * 1) * 256 +  __id[2]] = temp_2;
        
        outputs[(tx + bx * 1) * 256 +  __id[3]] = temp_3;
        
        outputs[(tx + bx * 1) * 256 +  __id[6]] = temp_6;
        
        outputs[(tx + bx * 1) * 256 +  __id[7]] = temp_7;
        
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