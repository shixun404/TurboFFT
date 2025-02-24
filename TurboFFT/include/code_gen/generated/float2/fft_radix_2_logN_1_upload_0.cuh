
#include "../../../TurboFFT_radix_2_template.h"
template<>
__global__ void fft_radix_2<float2, 1, 0, 0, 0, 0>(float2* inputs, float2* outputs, float2* twiddle, float2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    float2* shared = (float2*) ext_shared;
    int threadblock_per_SM = 256;
    int tb_gap = threadblock_per_SM * 108;
    int delta_bid = ((blockIdx.x / tb_gap) ==  (gridDim.x / tb_gap)) ? (gridDim.x % tb_gap) : tb_gap;
    float2 r[3];
    r[0].x = 1.0;
    r[0].y = 0.0;
    r[1].x = -0.5;
    r[1].y = -0.8660253882408142;
    r[2].x = -0.5;
    r[2].y = 0.8660253882408142;
    int j;
    int k;
    int global_j;
    int global_k;
    int data_id;
    int bs_id;
    int shared_offset_bs;
    int shared_offset_data;
    int bx;
    int tx;
    int offset;
    float2* gPtr;
    float2* shPtr;
    float2 rPtr[2];
    float2 rPtr_2[2];
    float2 rPtr_3[2];
    float2 rPtr_4[2];
    float2 tmp;
    float2 tmp_1;
    float2 tmp_2;
    float2 tmp_3;
    float2 angle;
    float2 delta_angle;
    j = 0;
    k = -1;
    global_j = 0;
    global_k = 0;
    data_id = 0;
    bs_id = 0;
    shared_offset_bs = 0;
    shared_offset_data = 0;
    bx = blockIdx.x;
    tx = threadIdx.x;
    offset = 0;
    gPtr = inputs;
    shPtr = shared;
    
    __syncthreads();
    int bid = 0;
    for(bid = (blockIdx.x / tb_gap) * tb_gap * thread_bs + blockIdx.x % tb_gap;
                bid_cnt < thread_bs && bid < (2 * BS + 64 - 1) / 64; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 32 * 1;
    
    gPtr += (bx % 1) * 2 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 1) * 32 * 2;
    bx = bx / 1;
    
    gPtr += tx % 32 * 2;
    
    gPtr += (bx % BS * 64);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        rPtr[1] = *(gPtr + 1);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[1]);
    turboFFT_ZSUB(rPtr[1], tmp, rPtr[1]);
    tmp = rPtr[1];
            
    bx = bid;
    tx = threadIdx.x;
    gPtr = outputs;
    
    gPtr += tx / 32 * 1;
    
    gPtr += (bx % 1) * 2 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 1) * 32 * 2;
    bx = bx / 1;
    
    gPtr += tx % 32 * 2;
    
    gPtr += (bx % BS * 64);
    
            *(gPtr + 0) = rPtr[0];
            
            *(gPtr + 1) = rPtr[1];
            
    }
    
}
