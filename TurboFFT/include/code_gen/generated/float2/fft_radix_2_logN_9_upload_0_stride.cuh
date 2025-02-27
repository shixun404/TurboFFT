
extern __shared__ float shared_mem[];
__global__ void fft_9_stride(float2* gPtr_1, float2* outputs, int threadblock_bs, int DY, int global_bs) {
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
    float2 rPtr[8];
    float2 rPtr_3[8];
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
    
    for(int bid_itr = 0; (bid_itr * gridDim.x + blockIdx.x) * threadblock_bs < global_bs; ++bid_itr){
    
gPtr = gPtr_1 + bid_itr * ((gridDim.x % (DY / threadblock_bs)) * threadblock_bs + (gridDim.x / (DY / threadblock_bs)) * DY * 512);
        
    bx = blockIdx.x;
    tx = threadIdx.x;
    
        gPtr += threadIdx.x % 64 * DY;
    
    gPtr += ((blockIdx.x % (DY / threadblock_bs)) * threadblock_bs + threadIdx.x / 64) + (blockIdx.x / (DY / threadblock_bs)) * DY * 512;
    

        rPtr[0] = *(gPtr + 0 * DY);
        
        rPtr[1] = *(gPtr + 64 * DY);
        
        rPtr[2] = *(gPtr + 128 * DY);
        
        rPtr[3] = *(gPtr + 192 * DY);
        
        rPtr[4] = *(gPtr + 256 * DY);
        
        rPtr[5] = *(gPtr + 320 * DY);
        
        rPtr[6] = *(gPtr + 384 * DY);
        
        rPtr[7] = *(gPtr + 448 * DY);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[5], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[6]);
    turboFFT_ZSUB(rPtr[6], tmp, rPtr[6]);
    tmp = rPtr[6];
    
    rPtr[6].y = -tmp.x;
    rPtr[6].x = tmp.y;
    
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
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
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[6]);
    turboFFT_ZSUB(rPtr[6], tmp, rPtr[6]);
    tmp = rPtr[6];
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
    rPtr[7].y = -tmp.x;
    rPtr[7].x = tmp.y;
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[1]);
    turboFFT_ZSUB(rPtr[1], tmp, rPtr[1]);
    tmp = rPtr[1];
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[3]);
    turboFFT_ZSUB(rPtr[3], tmp, rPtr[3]);
    tmp = rPtr[3];
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
    j = 0;
    offset  = 0;
    
    offset += ((threadIdx.x / 1) % 1) * 1;
    
    j = (threadIdx.x % 64) / 1;
    
    offset += ((threadIdx.x / 1) % 8) * 8;
    
    offset += ((threadIdx.x / 8) % 8) * 64;
    
    offset += (threadIdx.x / 64) * 512;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.012271846644580364f);
    delta_angle.y = __sinf(j * -0.012271846644580364f);
     
    angle.x = 1;
    angle.y = 0;
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
            rPtr_3[0] = rPtr[0];
    
            rPtr_3[1] = rPtr[4];
    
            rPtr_3[2] = rPtr[2];
    
            rPtr_3[3] = rPtr[6];
    
            rPtr_3[4] = rPtr[1];
    
            rPtr_3[5] = rPtr[5];
    
            rPtr_3[6] = rPtr[3];
    
            rPtr_3[7] = rPtr[7];
    
    shPtr[offset + 1 * ((0 + (threadIdx.x / 2)) % 8)] = rPtr_3[((0 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((0 + (threadIdx.x / 2)) % 8)] = rPtr_3[0];
     // shPtr[offset + 0] = rPtr[0];
    //  shPtr[offset + 0] = rPtr_3[0];
    
    shPtr[offset + 1 * ((1 + (threadIdx.x / 2)) % 8)] = rPtr_3[((1 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((1 + (threadIdx.x / 2)) % 8)] = rPtr_3[1];
     // shPtr[offset + 1] = rPtr[4];
    //  shPtr[offset + 1] = rPtr_3[1];
    
    shPtr[offset + 1 * ((2 + (threadIdx.x / 2)) % 8)] = rPtr_3[((2 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((2 + (threadIdx.x / 2)) % 8)] = rPtr_3[2];
     // shPtr[offset + 2] = rPtr[2];
    //  shPtr[offset + 2] = rPtr_3[2];
    
    shPtr[offset + 1 * ((3 + (threadIdx.x / 2)) % 8)] = rPtr_3[((3 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((3 + (threadIdx.x / 2)) % 8)] = rPtr_3[3];
     // shPtr[offset + 3] = rPtr[6];
    //  shPtr[offset + 3] = rPtr_3[3];
    
    shPtr[offset + 1 * ((4 + (threadIdx.x / 2)) % 8)] = rPtr_3[((4 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((4 + (threadIdx.x / 2)) % 8)] = rPtr_3[4];
     // shPtr[offset + 4] = rPtr[1];
    //  shPtr[offset + 4] = rPtr_3[4];
    
    shPtr[offset + 1 * ((5 + (threadIdx.x / 2)) % 8)] = rPtr_3[((5 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((5 + (threadIdx.x / 2)) % 8)] = rPtr_3[5];
     // shPtr[offset + 5] = rPtr[5];
    //  shPtr[offset + 5] = rPtr_3[5];
    
    shPtr[offset + 1 * ((6 + (threadIdx.x / 2)) % 8)] = rPtr_3[((6 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((6 + (threadIdx.x / 2)) % 8)] = rPtr_3[6];
     // shPtr[offset + 6] = rPtr[3];
    //  shPtr[offset + 6] = rPtr_3[6];
    
    shPtr[offset + 1 * ((7 + (threadIdx.x / 2)) % 8)] = rPtr_3[((7 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 1 * ((7 + (threadIdx.x / 2)) % 8)] = rPtr_3[7];
     // shPtr[offset + 7] = rPtr[7];
    //  shPtr[offset + 7] = rPtr_3[7];
    
    offset = 0;
    offset += tx % 64 + tx / 64 * 512;
    
    __syncthreads();
    
    rPtr[0] = shPtr[offset + 0];
    
    rPtr[1] = shPtr[offset + 64];
    
    rPtr[2] = shPtr[offset + 128];
    
    rPtr[3] = shPtr[offset + 192];
    
    rPtr[4] = shPtr[offset + 256];
    
    rPtr[5] = shPtr[offset + 320];
    
    rPtr[6] = shPtr[offset + 384];
    
    rPtr[7] = shPtr[offset + 448];
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[5], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[6]);
    turboFFT_ZSUB(rPtr[6], tmp, rPtr[6]);
    tmp = rPtr[6];
    
    rPtr[6].y = -tmp.x;
    rPtr[6].x = tmp.y;
    
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
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
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[6]);
    turboFFT_ZSUB(rPtr[6], tmp, rPtr[6]);
    tmp = rPtr[6];
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
    rPtr[7].y = -tmp.x;
    rPtr[7].x = tmp.y;
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[1]);
    turboFFT_ZSUB(rPtr[1], tmp, rPtr[1]);
    tmp = rPtr[1];
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[3]);
    turboFFT_ZSUB(rPtr[3], tmp, rPtr[3]);
    tmp = rPtr[3];
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
    j = 0;
    offset  = 0;
    
    offset += ((threadIdx.x / 1) % 1) * 1;
    
    offset += ((threadIdx.x / 1) % 8) * 1;
    
    j = (threadIdx.x % 64) / 8;
    
    offset += ((threadIdx.x / 8) % 8) * 64;
    
    offset += (threadIdx.x / 64) * 512;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.09817477315664291f);
    delta_angle.y = __sinf(j * -0.09817477315664291f);
     
    angle.x = 1;
    angle.y = 0;
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
            rPtr_3[0] = rPtr[0];
    
            rPtr_3[1] = rPtr[4];
    
            rPtr_3[2] = rPtr[2];
    
            rPtr_3[3] = rPtr[6];
    
            rPtr_3[4] = rPtr[1];
    
            rPtr_3[5] = rPtr[5];
    
            rPtr_3[6] = rPtr[3];
    
            rPtr_3[7] = rPtr[7];
    
    shPtr[offset + 8 * ((0 + (threadIdx.x / 2)) % 8)] = rPtr_3[((0 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((0 + (threadIdx.x / 2)) % 8)] = rPtr_3[0];
     // shPtr[offset + 0] = rPtr[0];
     // shPtr[offset + 0] = rPtr_3[0];
    
    shPtr[offset + 8 * ((1 + (threadIdx.x / 2)) % 8)] = rPtr_3[((1 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((1 + (threadIdx.x / 2)) % 8)] = rPtr_3[1];
     // shPtr[offset + 8] = rPtr[4];
     // shPtr[offset + 8] = rPtr_3[1];
    
    shPtr[offset + 8 * ((2 + (threadIdx.x / 2)) % 8)] = rPtr_3[((2 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((2 + (threadIdx.x / 2)) % 8)] = rPtr_3[2];
     // shPtr[offset + 16] = rPtr[2];
     // shPtr[offset + 16] = rPtr_3[2];
    
    shPtr[offset + 8 * ((3 + (threadIdx.x / 2)) % 8)] = rPtr_3[((3 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((3 + (threadIdx.x / 2)) % 8)] = rPtr_3[3];
     // shPtr[offset + 24] = rPtr[6];
     // shPtr[offset + 24] = rPtr_3[3];
    
    shPtr[offset + 8 * ((4 + (threadIdx.x / 2)) % 8)] = rPtr_3[((4 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((4 + (threadIdx.x / 2)) % 8)] = rPtr_3[4];
     // shPtr[offset + 32] = rPtr[1];
     // shPtr[offset + 32] = rPtr_3[4];
    
    shPtr[offset + 8 * ((5 + (threadIdx.x / 2)) % 8)] = rPtr_3[((5 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((5 + (threadIdx.x / 2)) % 8)] = rPtr_3[5];
     // shPtr[offset + 40] = rPtr[5];
     // shPtr[offset + 40] = rPtr_3[5];
    
    shPtr[offset + 8 * ((6 + (threadIdx.x / 2)) % 8)] = rPtr_3[((6 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((6 + (threadIdx.x / 2)) % 8)] = rPtr_3[6];
     // shPtr[offset + 48] = rPtr[3];
     // shPtr[offset + 48] = rPtr_3[6];
    
    shPtr[offset + 8 * ((7 + (threadIdx.x / 2)) % 8)] = rPtr_3[((7 + (threadIdx.x / 2)) % 8)];
    // shPtr[offset + 8 * ((7 + (threadIdx.x / 2)) % 8)] = rPtr_3[7];
     // shPtr[offset + 56] = rPtr[7];
     // shPtr[offset + 56] = rPtr_3[7];
    
    offset = 0;
    offset += tx % 64 + tx / 64 * 512;
    
    __syncthreads();
    
    rPtr[0] = shPtr[offset + 0];
    
    rPtr[1] = shPtr[offset + 64];
    
    rPtr[2] = shPtr[offset + 128];
    
    rPtr[3] = shPtr[offset + 192];
    
    rPtr[4] = shPtr[offset + 256];
    
    rPtr[5] = shPtr[offset + 320];
    
    rPtr[6] = shPtr[offset + 384];
    
    rPtr[7] = shPtr[offset + 448];
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[5], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[6]);
    turboFFT_ZSUB(rPtr[6], tmp, rPtr[6]);
    tmp = rPtr[6];
    
    rPtr[6].y = -tmp.x;
    rPtr[6].x = tmp.y;
    
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
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
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[6]);
    turboFFT_ZSUB(rPtr[6], tmp, rPtr[6]);
    tmp = rPtr[6];
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
    
    rPtr[7].y = -tmp.x;
    rPtr[7].x = tmp.y;
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[1]);
    turboFFT_ZSUB(rPtr[1], tmp, rPtr[1]);
    tmp = rPtr[1];
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[3]);
    turboFFT_ZSUB(rPtr[3], tmp, rPtr[3]);
    tmp = rPtr[3];
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[7]);
    turboFFT_ZSUB(rPtr[7], tmp, rPtr[7]);
    tmp = rPtr[7];
            
    bx = blockIdx.x;
    tx = threadIdx.x;
    gPtr = outputs + bid_itr * ((gridDim.x % (DY / threadblock_bs)) * threadblock_bs + (gridDim.x / (DY / threadblock_bs)) * DY * 64);
            gPtr += threadIdx.x % 64 * DY;
    
    gPtr += ((blockIdx.x % (DY / threadblock_bs)) * threadblock_bs + threadIdx.x / 64) + (blockIdx.x / (DY / threadblock_bs)) * DY * 64;
    
                rPtr_3[0] = rPtr[0];
        
                rPtr_3[1] = rPtr[4];
        
                rPtr_3[2] = rPtr[2];
        
                rPtr_3[3] = rPtr[6];
        
                rPtr_3[4] = rPtr[1];
        
                rPtr_3[5] = rPtr[5];
        
                rPtr_3[6] = rPtr[3];
        
                rPtr_3[7] = rPtr[7];
        
                    #pragma unroll
                    for(int i = 0; i < (1); ++i)
                    *(gPtr + i * 64 * DY) = rPtr_3[i];
            
}
}
