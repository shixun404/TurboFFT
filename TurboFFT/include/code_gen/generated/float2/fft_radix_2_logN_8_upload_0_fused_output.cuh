
extern __shared__ float shared_mem[];
__device__ void fft_8_fused_output(float2* inputs, float2* outputs, float2* sFFT, int stride) {
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
    float2 rPtr[16];
    float2 rPtr_3[16];
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
    gPtr = inputs;
    shPtr = sFFT;
    
    int bid = 0;
            
    bx = blockIdx.x;
    tx = threadIdx.x;
    
        if(threadIdx.x < THREADBLOCK_K * 16 )
        {
        gPtr = inputs;
        gPtr += threadIdx.x % 16;
        gPtr += (threadIdx.x / 16) * 256;


        rPtr[0] = {0, 0};
        rPtr[1] = {0, 0};
        rPtr[2] = {0, 0};
        rPtr[3] = {0, 0};
        rPtr[4] = {0, 0};
        rPtr[5] = {0, 0};
        rPtr[6] = {0, 0};
        rPtr[7] = {0, 0};
        rPtr[8] = {0, 0};
        rPtr[9] = {0, 0};
        rPtr[10] = {0, 0};
        rPtr[11] = {0, 0};
        rPtr[12] = {0, 0};
        rPtr[13] = {0, 0};
        rPtr[14] = {0, 0};
        rPtr[15] = {0, 0};
                #pragma unroll
                for(int i = 0; i < (THREADBLOCK_M / 16); ++i)
                rPtr[i] = *(gPtr + i * 16);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867f;
        angle.y = -0.3826834323650898f;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[11], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    rPtr[12].y = -tmp.x;
    rPtr[12].x = tmp.y;
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = -0.3826834323650897f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867f;
        angle.y = -0.3826834323650899f;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
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
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
    rPtr[14].y = -tmp.x;
    rPtr[14].x = tmp.y;
    
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
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
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
    rPtr[11].y = -tmp.x;
    rPtr[11].x = tmp.y;
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
    rPtr[15].y = -tmp.x;
    rPtr[15].x = tmp.y;
    
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
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
    j = 0;
    offset  = 0;
    
    offset += ((threadIdx.x / 1) % 1) * 1;
    
    j = (threadIdx.x % 16) / 1;
    
    offset += ((threadIdx.x / 1) % 16) * 16;
    
    offset += (threadIdx.x / 16) * 256;
    
    }
    __syncthreads();
    
    if(threadIdx.x < THREADBLOCK_K * 16 )
    {
    delta_angle.x = __cosf(j * -0.02454369328916073f);
    delta_angle.y = __sinf(j * -0.02454369328916073f);
     
    angle.x = 1;
    angle.y = 0;
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
            rPtr_3[0] = rPtr[0];
    
            rPtr_3[1] = rPtr[8];
    
            rPtr_3[2] = rPtr[4];
    
            rPtr_3[3] = rPtr[12];
    
            rPtr_3[4] = rPtr[2];
    
            rPtr_3[5] = rPtr[10];
    
            rPtr_3[6] = rPtr[6];
    
            rPtr_3[7] = rPtr[14];
    
            rPtr_3[8] = rPtr[1];
    
            rPtr_3[9] = rPtr[9];
    
            rPtr_3[10] = rPtr[5];
    
            rPtr_3[11] = rPtr[13];
    
            rPtr_3[12] = rPtr[3];
    
            rPtr_3[13] = rPtr[11];
    
            rPtr_3[14] = rPtr[7];
    
            rPtr_3[15] = rPtr[15];
    
    shPtr[offset + 1 * ((0 + (threadIdx.x / 1)) % 16)] = rPtr_3[((0 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((0 + (threadIdx.x / 1)) % 16)] = rPtr_3[0];
     // shPtr[offset + 0] = rPtr[0];
    //  shPtr[offset + 0] = rPtr_3[0];
    
    shPtr[offset + 1 * ((1 + (threadIdx.x / 1)) % 16)] = rPtr_3[((1 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((1 + (threadIdx.x / 1)) % 16)] = rPtr_3[1];
     // shPtr[offset + 1] = rPtr[8];
    //  shPtr[offset + 1] = rPtr_3[1];
    
    shPtr[offset + 1 * ((2 + (threadIdx.x / 1)) % 16)] = rPtr_3[((2 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((2 + (threadIdx.x / 1)) % 16)] = rPtr_3[2];
     // shPtr[offset + 2] = rPtr[4];
    //  shPtr[offset + 2] = rPtr_3[2];
    
    shPtr[offset + 1 * ((3 + (threadIdx.x / 1)) % 16)] = rPtr_3[((3 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((3 + (threadIdx.x / 1)) % 16)] = rPtr_3[3];
     // shPtr[offset + 3] = rPtr[12];
    //  shPtr[offset + 3] = rPtr_3[3];
    
    shPtr[offset + 1 * ((4 + (threadIdx.x / 1)) % 16)] = rPtr_3[((4 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((4 + (threadIdx.x / 1)) % 16)] = rPtr_3[4];
     // shPtr[offset + 4] = rPtr[2];
    //  shPtr[offset + 4] = rPtr_3[4];
    
    shPtr[offset + 1 * ((5 + (threadIdx.x / 1)) % 16)] = rPtr_3[((5 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((5 + (threadIdx.x / 1)) % 16)] = rPtr_3[5];
     // shPtr[offset + 5] = rPtr[10];
    //  shPtr[offset + 5] = rPtr_3[5];
    
    shPtr[offset + 1 * ((6 + (threadIdx.x / 1)) % 16)] = rPtr_3[((6 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((6 + (threadIdx.x / 1)) % 16)] = rPtr_3[6];
     // shPtr[offset + 6] = rPtr[6];
    //  shPtr[offset + 6] = rPtr_3[6];
    
    shPtr[offset + 1 * ((7 + (threadIdx.x / 1)) % 16)] = rPtr_3[((7 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((7 + (threadIdx.x / 1)) % 16)] = rPtr_3[7];
     // shPtr[offset + 7] = rPtr[14];
    //  shPtr[offset + 7] = rPtr_3[7];
    
    shPtr[offset + 1 * ((8 + (threadIdx.x / 1)) % 16)] = rPtr_3[((8 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((8 + (threadIdx.x / 1)) % 16)] = rPtr_3[8];
     // shPtr[offset + 8] = rPtr[1];
    //  shPtr[offset + 8] = rPtr_3[8];
    
    shPtr[offset + 1 * ((9 + (threadIdx.x / 1)) % 16)] = rPtr_3[((9 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((9 + (threadIdx.x / 1)) % 16)] = rPtr_3[9];
     // shPtr[offset + 9] = rPtr[9];
    //  shPtr[offset + 9] = rPtr_3[9];
    
    shPtr[offset + 1 * ((10 + (threadIdx.x / 1)) % 16)] = rPtr_3[((10 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((10 + (threadIdx.x / 1)) % 16)] = rPtr_3[10];
     // shPtr[offset + 10] = rPtr[5];
    //  shPtr[offset + 10] = rPtr_3[10];
    
    shPtr[offset + 1 * ((11 + (threadIdx.x / 1)) % 16)] = rPtr_3[((11 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((11 + (threadIdx.x / 1)) % 16)] = rPtr_3[11];
     // shPtr[offset + 11] = rPtr[13];
    //  shPtr[offset + 11] = rPtr_3[11];
    
    shPtr[offset + 1 * ((12 + (threadIdx.x / 1)) % 16)] = rPtr_3[((12 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((12 + (threadIdx.x / 1)) % 16)] = rPtr_3[12];
     // shPtr[offset + 12] = rPtr[3];
    //  shPtr[offset + 12] = rPtr_3[12];
    
    shPtr[offset + 1 * ((13 + (threadIdx.x / 1)) % 16)] = rPtr_3[((13 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((13 + (threadIdx.x / 1)) % 16)] = rPtr_3[13];
     // shPtr[offset + 13] = rPtr[11];
    //  shPtr[offset + 13] = rPtr_3[13];
    
    shPtr[offset + 1 * ((14 + (threadIdx.x / 1)) % 16)] = rPtr_3[((14 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((14 + (threadIdx.x / 1)) % 16)] = rPtr_3[14];
     // shPtr[offset + 14] = rPtr[7];
    //  shPtr[offset + 14] = rPtr_3[14];
    
    shPtr[offset + 1 * ((15 + (threadIdx.x / 1)) % 16)] = rPtr_3[((15 + (threadIdx.x / 1)) % 16)];
    // shPtr[offset + 1 * ((15 + (threadIdx.x / 1)) % 16)] = rPtr_3[15];
     // shPtr[offset + 15] = rPtr[15];
    //  shPtr[offset + 15] = rPtr_3[15];
    
    offset = 0;
    offset += tx % 16 + tx / 16 * 256;
    
    }
    __syncthreads();    
    
    if(threadIdx.x < THREADBLOCK_K * 16 )
    {
    rPtr[0] = shPtr[offset + 0];
    
    rPtr[1] = shPtr[offset + 16];
    
    rPtr[2] = shPtr[offset + 32];
    
    rPtr[3] = shPtr[offset + 48];
    
    rPtr[4] = shPtr[offset + 64];
    
    rPtr[5] = shPtr[offset + 80];
    
    rPtr[6] = shPtr[offset + 96];
    
    rPtr[7] = shPtr[offset + 112];
    
    rPtr[8] = shPtr[offset + 128];
    
    rPtr[9] = shPtr[offset + 144];
    
    rPtr[10] = shPtr[offset + 160];
    
    rPtr[11] = shPtr[offset + 176];
    
    rPtr[12] = shPtr[offset + 192];
    
    rPtr[13] = shPtr[offset + 208];
    
    rPtr[14] = shPtr[offset + 224];
    
    rPtr[15] = shPtr[offset + 240];
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867f;
        angle.y = -0.3826834323650898f;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[11], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    rPtr[12].y = -tmp.x;
    rPtr[12].x = tmp.y;
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = -0.3826834323650897f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867f;
        angle.y = -0.3826834323650899f;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
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
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
    rPtr[14].y = -tmp.x;
    rPtr[14].x = tmp.y;
    
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
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
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
    rPtr[11].y = -tmp.x;
    rPtr[11].x = tmp.y;
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
    rPtr[15].y = -tmp.x;
    rPtr[15].x = tmp.y;
    
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
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
            
    bx = blockIdx.x;
    tx = threadIdx.x;
    gPtr = outputs;
        gPtr += threadIdx.x % 16;
    
        gPtr += (threadIdx.x / 16) * stride;
    
                    *(gPtr + 0) = rPtr[0];
                    
                    *(gPtr + 16) = rPtr[8];
                    
                    *(gPtr + 32) = rPtr[4];
                    
                    *(gPtr + 48) = rPtr[12];
                    
                    *(gPtr + 64) = rPtr[2];
                    
                    *(gPtr + 80) = rPtr[10];
                    
                    *(gPtr + 96) = rPtr[6];
                    
                    *(gPtr + 112) = rPtr[14];
                    
                    *(gPtr + 128) = rPtr[1];
                    
                    *(gPtr + 144) = rPtr[9];
                    
                    *(gPtr + 160) = rPtr[5];
                    
                    *(gPtr + 176) = rPtr[13];
                    
                    *(gPtr + 192) = rPtr[3];
                    
                    *(gPtr + 208) = rPtr[11];
                    
                    *(gPtr + 224) = rPtr[7];
                    
                    *(gPtr + 240) = rPtr[15];
                    }
}
