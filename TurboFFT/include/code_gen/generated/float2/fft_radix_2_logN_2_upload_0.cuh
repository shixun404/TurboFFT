
extern __shared__ float shared_mem[];
__global__ void fft_2(float2* gPtr_1, float2* outputs, int threadblock_bs) {
    int bid_cnt = 0;
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
    float2 rPtr[4];
    float2 rPtr_3[4];
    float2 tmp;
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
    gPtr = gPtr_1;
    shPtr = (float2*) shared_mem;
    
    int bid = 0;
            
    bx = blockIdx.x;
    tx = threadIdx.x;
    
        gPtr += threadIdx.x % 1;
    
    gPtr += (blockIdx.x * threadblock_bs + threadIdx.x / 1) * 4;
    

        rPtr[0] = *(gPtr + 0);
        
        rPtr[1] = *(gPtr + 1);
        
        rPtr[2] = *(gPtr + 2);
        
        rPtr[3] = *(gPtr + 3);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[2]);
    turboFFT_ZSUB(rPtr[2], tmp, rPtr[2]);
    tmp = rPtr[2];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[3]);
    turboFFT_ZSUB(rPtr[3], tmp, rPtr[3]);
    tmp = rPtr[3];
    
    rPtr[3].y = -tmp.x;
    rPtr[3].x = tmp.y;
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[1]);
    turboFFT_ZSUB(rPtr[1], tmp, rPtr[1]);
    tmp = rPtr[1];
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[3]);
    turboFFT_ZSUB(rPtr[3], tmp, rPtr[3]);
    tmp = rPtr[3];
            
    bx = blockIdx.x;
    tx = threadIdx.x;
    gPtr = outputs;
    
        gPtr += threadIdx.x % 1;
    
    gPtr += (blockIdx.x * threadblock_bs + threadIdx.x / 1) * 4;
    

            *(gPtr + 0) = rPtr[0];
            
            *(gPtr + 1) = rPtr[2];
            
            *(gPtr + 2) = rPtr[1];
            
            *(gPtr + 3) = rPtr[3];
            
}
