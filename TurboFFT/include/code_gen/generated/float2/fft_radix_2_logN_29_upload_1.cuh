
#include "../../../TurboFFT_radix_2_template.h"
template<>
__global__ void fft_radix_2<float2, 29, 1, 0, 0, 0>(float2* inputs, float2* outputs, float2* twiddle, float2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    float2* shared = (float2*) ext_shared;
    int threadblock_per_SM = 2;
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
    float2 rPtr[32];
    float2 rPtr_2[32];
    float2 rPtr_3[32];
    float2 rPtr_4[32];
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
                bid_cnt < thread_bs && bid < (536870912 * BS + 8192 - 1) / 8192; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += (bx % 128) * 8 * 1;
    bx = bx / 128;
    
    gPtr += tx % 8 * 1;
    
    gPtr += tx / 8 * 1024;
    
    gPtr += (bx % 1) * 1024 * 1024;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 1048576;
    bx = bx / 512;
    
    gPtr += (bx % BS * 536870912);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        rPtr[1] = *(gPtr + 32768);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
        rPtr[2] = *(gPtr + 65536);
        rPtr_3[2].x += rPtr[2].x;
        rPtr_3[2].y += rPtr[2].y;
        
        rPtr[3] = *(gPtr + 98304);
        rPtr_3[3].x += rPtr[3].x;
        rPtr_3[3].y += rPtr[3].y;
        
        rPtr[4] = *(gPtr + 131072);
        rPtr_3[4].x += rPtr[4].x;
        rPtr_3[4].y += rPtr[4].y;
        
        rPtr[5] = *(gPtr + 163840);
        rPtr_3[5].x += rPtr[5].x;
        rPtr_3[5].y += rPtr[5].y;
        
        rPtr[6] = *(gPtr + 196608);
        rPtr_3[6].x += rPtr[6].x;
        rPtr_3[6].y += rPtr[6].y;
        
        rPtr[7] = *(gPtr + 229376);
        rPtr_3[7].x += rPtr[7].x;
        rPtr_3[7].y += rPtr[7].y;
        
        rPtr[8] = *(gPtr + 262144);
        rPtr_3[8].x += rPtr[8].x;
        rPtr_3[8].y += rPtr[8].y;
        
        rPtr[9] = *(gPtr + 294912);
        rPtr_3[9].x += rPtr[9].x;
        rPtr_3[9].y += rPtr[9].y;
        
        rPtr[10] = *(gPtr + 327680);
        rPtr_3[10].x += rPtr[10].x;
        rPtr_3[10].y += rPtr[10].y;
        
        rPtr[11] = *(gPtr + 360448);
        rPtr_3[11].x += rPtr[11].x;
        rPtr_3[11].y += rPtr[11].y;
        
        rPtr[12] = *(gPtr + 393216);
        rPtr_3[12].x += rPtr[12].x;
        rPtr_3[12].y += rPtr[12].y;
        
        rPtr[13] = *(gPtr + 425984);
        rPtr_3[13].x += rPtr[13].x;
        rPtr_3[13].y += rPtr[13].y;
        
        rPtr[14] = *(gPtr + 458752);
        rPtr_3[14].x += rPtr[14].x;
        rPtr_3[14].y += rPtr[14].y;
        
        rPtr[15] = *(gPtr + 491520);
        rPtr_3[15].x += rPtr[15].x;
        rPtr_3[15].y += rPtr[15].y;
        
        rPtr[16] = *(gPtr + 524288);
        rPtr_3[16].x += rPtr[16].x;
        rPtr_3[16].y += rPtr[16].y;
        
        rPtr[17] = *(gPtr + 557056);
        rPtr_3[17].x += rPtr[17].x;
        rPtr_3[17].y += rPtr[17].y;
        
        rPtr[18] = *(gPtr + 589824);
        rPtr_3[18].x += rPtr[18].x;
        rPtr_3[18].y += rPtr[18].y;
        
        rPtr[19] = *(gPtr + 622592);
        rPtr_3[19].x += rPtr[19].x;
        rPtr_3[19].y += rPtr[19].y;
        
        rPtr[20] = *(gPtr + 655360);
        rPtr_3[20].x += rPtr[20].x;
        rPtr_3[20].y += rPtr[20].y;
        
        rPtr[21] = *(gPtr + 688128);
        rPtr_3[21].x += rPtr[21].x;
        rPtr_3[21].y += rPtr[21].y;
        
        rPtr[22] = *(gPtr + 720896);
        rPtr_3[22].x += rPtr[22].x;
        rPtr_3[22].y += rPtr[22].y;
        
        rPtr[23] = *(gPtr + 753664);
        rPtr_3[23].x += rPtr[23].x;
        rPtr_3[23].y += rPtr[23].y;
        
        rPtr[24] = *(gPtr + 786432);
        rPtr_3[24].x += rPtr[24].x;
        rPtr_3[24].y += rPtr[24].y;
        
        rPtr[25] = *(gPtr + 819200);
        rPtr_3[25].x += rPtr[25].x;
        rPtr_3[25].y += rPtr[25].y;
        
        rPtr[26] = *(gPtr + 851968);
        rPtr_3[26].x += rPtr[26].x;
        rPtr_3[26].y += rPtr[26].y;
        
        rPtr[27] = *(gPtr + 884736);
        rPtr_3[27].x += rPtr[27].x;
        rPtr_3[27].y += rPtr[27].y;
        
        rPtr[28] = *(gPtr + 917504);
        rPtr_3[28].x += rPtr[28].x;
        rPtr_3[28].y += rPtr[28].y;
        
        rPtr[29] = *(gPtr + 950272);
        rPtr_3[29].x += rPtr[29].x;
        rPtr_3[29].y += rPtr[29].y;
        
        rPtr[30] = *(gPtr + 983040);
        rPtr_3[30].x += rPtr[30].x;
        rPtr_3[30].y += rPtr[30].y;
        
        rPtr[31] = *(gPtr + 1015808);
        rPtr_3[31].x += rPtr[31].x;
        rPtr_3[31].y += rPtr[31].y;
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[16]);
    turboFFT_ZSUB(rPtr[16], tmp, rPtr[16]);
    tmp = rPtr[16];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[17]);
    turboFFT_ZSUB(rPtr[17], tmp, rPtr[17]);
    tmp = rPtr[17];
    
        angle.x = 0.9807852804032304f;
        angle.y = -0.19509032201612825f;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867f;
        angle.y = -0.3826834323650898f;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452f;
        angle.y = -0.5555702330196022f;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023f;
        angle.y = -0.8314696123025452f;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833f;
        angle.y = -0.9807852804032304f;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    rPtr[24].y = -tmp.x;
    rPtr[24].x = tmp.y;
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = -0.1950903220161282f;
        angle.y = -0.9807852804032304f;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602f;
        angle.y = -0.8314696123025455f;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453f;
        angle.y = -0.5555702330196022f;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867f;
        angle.y = -0.3826834323650899f;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304f;
        angle.y = -0.1950903220161286f;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
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
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867f;
        angle.y = -0.3826834323650898f;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[20];
    turboFFT_ZADD(rPtr[20], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    rPtr[28].y = -tmp.x;
    rPtr[28].x = tmp.y;
    
    tmp = rPtr[21];
    turboFFT_ZADD(rPtr[21], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.3826834323650897f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867f;
        angle.y = -0.3826834323650899f;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
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
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    rPtr[22].y = -tmp.x;
    rPtr[22].x = tmp.y;
    
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[26];
    turboFFT_ZADD(rPtr[26], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    rPtr[30].y = -tmp.x;
    rPtr[30].x = tmp.y;
    
    tmp = rPtr[27];
    turboFFT_ZADD(rPtr[27], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
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
    
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    rPtr[19].y = -tmp.x;
    rPtr[19].x = tmp.y;
    
    tmp = rPtr[20];
    turboFFT_ZADD(rPtr[20], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    tmp = rPtr[21];
    turboFFT_ZADD(rPtr[21], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    rPtr[23].y = -tmp.x;
    rPtr[23].x = tmp.y;
    
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    rPtr[27].y = -tmp.x;
    rPtr[27].x = tmp.y;
    
    tmp = rPtr[28];
    turboFFT_ZADD(rPtr[28], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    tmp = rPtr[29];
    turboFFT_ZADD(rPtr[29], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
    rPtr[31].y = -tmp.x;
    rPtr[31].x = tmp.y;
    
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
    
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[17]);
    turboFFT_ZSUB(rPtr[17], tmp, rPtr[17]);
    tmp = rPtr[17];
    
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    tmp = rPtr[20];
    turboFFT_ZADD(rPtr[20], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
    tmp = rPtr[26];
    turboFFT_ZADD(rPtr[26], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    tmp = rPtr[28];
    turboFFT_ZADD(rPtr[28], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
    tmp = rPtr[30];
    turboFFT_ZADD(rPtr[30], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
    j = 0;
    offset  = 0;
    
    offset += ((tx / 1) % 8) * 1;
    
    j = tx / 8;
    
    offset += ((tx / 8) % 32) * 256;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.006135923322290182f);
    delta_angle.y = __sinf(j * -0.006135923322290182f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 8] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 16] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 24] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 32] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 40] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 48] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 56] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 64] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 72] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 80] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 88] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 96] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 104] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 112] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 120] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 128] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 136] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 144] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 152] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 160] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 168] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 176] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 184] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 192] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 200] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 208] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 216] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 224] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 232] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 240] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 248] = rPtr[31];
    
    offset = 0;
    offset += tx;
    
    __syncthreads();
    
    rPtr[0] = shPtr[offset + 0];
    
    rPtr[1] = shPtr[offset + 256];
    
    rPtr[2] = shPtr[offset + 512];
    
    rPtr[3] = shPtr[offset + 768];
    
    rPtr[4] = shPtr[offset + 1024];
    
    rPtr[5] = shPtr[offset + 1280];
    
    rPtr[6] = shPtr[offset + 1536];
    
    rPtr[7] = shPtr[offset + 1792];
    
    rPtr[8] = shPtr[offset + 2048];
    
    rPtr[9] = shPtr[offset + 2304];
    
    rPtr[10] = shPtr[offset + 2560];
    
    rPtr[11] = shPtr[offset + 2816];
    
    rPtr[12] = shPtr[offset + 3072];
    
    rPtr[13] = shPtr[offset + 3328];
    
    rPtr[14] = shPtr[offset + 3584];
    
    rPtr[15] = shPtr[offset + 3840];
    
    rPtr[16] = shPtr[offset + 4096];
    
    rPtr[17] = shPtr[offset + 4352];
    
    rPtr[18] = shPtr[offset + 4608];
    
    rPtr[19] = shPtr[offset + 4864];
    
    rPtr[20] = shPtr[offset + 5120];
    
    rPtr[21] = shPtr[offset + 5376];
    
    rPtr[22] = shPtr[offset + 5632];
    
    rPtr[23] = shPtr[offset + 5888];
    
    rPtr[24] = shPtr[offset + 6144];
    
    rPtr[25] = shPtr[offset + 6400];
    
    rPtr[26] = shPtr[offset + 6656];
    
    rPtr[27] = shPtr[offset + 6912];
    
    rPtr[28] = shPtr[offset + 7168];
    
    rPtr[29] = shPtr[offset + 7424];
    
    rPtr[30] = shPtr[offset + 7680];
    
    rPtr[31] = shPtr[offset + 7936];
    
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[16]);
    turboFFT_ZSUB(rPtr[16], tmp, rPtr[16]);
    tmp = rPtr[16];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[17]);
    turboFFT_ZSUB(rPtr[17], tmp, rPtr[17]);
    tmp = rPtr[17];
    
        angle.x = 0.9807852804032304f;
        angle.y = -0.19509032201612825f;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867f;
        angle.y = -0.3826834323650898f;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452f;
        angle.y = -0.5555702330196022f;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023f;
        angle.y = -0.8314696123025452f;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833f;
        angle.y = -0.9807852804032304f;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    rPtr[24].y = -tmp.x;
    rPtr[24].x = tmp.y;
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = -0.1950903220161282f;
        angle.y = -0.9807852804032304f;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602f;
        angle.y = -0.8314696123025455f;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453f;
        angle.y = -0.5555702330196022f;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867f;
        angle.y = -0.3826834323650899f;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304f;
        angle.y = -0.1950903220161286f;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
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
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867f;
        angle.y = -0.3826834323650898f;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[20];
    turboFFT_ZADD(rPtr[20], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    rPtr[28].y = -tmp.x;
    rPtr[28].x = tmp.y;
    
    tmp = rPtr[21];
    turboFFT_ZADD(rPtr[21], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.3826834323650897f;
        angle.y = -0.9238795325112867f;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867f;
        angle.y = -0.3826834323650899f;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
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
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    rPtr[22].y = -tmp.x;
    rPtr[22].x = tmp.y;
    
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476f;
        angle.y = -0.7071067811865475f;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[26];
    turboFFT_ZADD(rPtr[26], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    rPtr[30].y = -tmp.x;
    rPtr[30].x = tmp.y;
    
    tmp = rPtr[27];
    turboFFT_ZADD(rPtr[27], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.7071067811865475f;
        angle.y = -0.7071067811865476f;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
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
    
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    rPtr[19].y = -tmp.x;
    rPtr[19].x = tmp.y;
    
    tmp = rPtr[20];
    turboFFT_ZADD(rPtr[20], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    tmp = rPtr[21];
    turboFFT_ZADD(rPtr[21], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    rPtr[23].y = -tmp.x;
    rPtr[23].x = tmp.y;
    
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    rPtr[27].y = -tmp.x;
    rPtr[27].x = tmp.y;
    
    tmp = rPtr[28];
    turboFFT_ZADD(rPtr[28], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    tmp = rPtr[29];
    turboFFT_ZADD(rPtr[29], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
    rPtr[31].y = -tmp.x;
    rPtr[31].x = tmp.y;
    
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
    
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[17]);
    turboFFT_ZSUB(rPtr[17], tmp, rPtr[17]);
    tmp = rPtr[17];
    
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    tmp = rPtr[20];
    turboFFT_ZADD(rPtr[20], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
    tmp = rPtr[26];
    turboFFT_ZADD(rPtr[26], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    tmp = rPtr[28];
    turboFFT_ZADD(rPtr[28], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
    tmp = rPtr[30];
    turboFFT_ZADD(rPtr[30], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
            
    bx = bid;
    tx = threadIdx.x;
    gPtr = outputs;
    global_j = 0;
    global_k = 0;
    
    global_j += (bx % 128) * 8 * 1;
    
    global_j += (tx % 8) * 1;
    
    gPtr += (bx % 128) * 8 * 1;
    bx = bx / 128;
    
    gPtr += tx % 8 * 1;
    
    gPtr += tx / 8 * 1024;
    
    gPtr += (bx % 1) * 1024 * 1024;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 1048576;
    bx = bx / 512;
    
    gPtr += (bx % BS * 536870912);
    
    global_k += tx / 8;
    
        delta_angle.x = __cosf(global_j *  -0.0001917476038215682f);
        delta_angle.y = __sinf(global_j *  -0.0001917476038215682f);
        angle.x = __cosf( global_j * global_k * -5.992112452678286e-06f);
        angle.y = __sinf( global_j * global_k * -5.992112452678286e-06f);
        
            tmp = rPtr[0];
            turboFFT_ZMUL(rPtr[0], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[16];
            turboFFT_ZMUL(rPtr[16], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[8];
            turboFFT_ZMUL(rPtr[8], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[24];
            turboFFT_ZMUL(rPtr[24], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[4];
            turboFFT_ZMUL(rPtr[4], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[20];
            turboFFT_ZMUL(rPtr[20], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[12];
            turboFFT_ZMUL(rPtr[12], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[28];
            turboFFT_ZMUL(rPtr[28], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[2];
            turboFFT_ZMUL(rPtr[2], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[18];
            turboFFT_ZMUL(rPtr[18], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[10];
            turboFFT_ZMUL(rPtr[10], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[26];
            turboFFT_ZMUL(rPtr[26], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[6];
            turboFFT_ZMUL(rPtr[6], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[22];
            turboFFT_ZMUL(rPtr[22], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[14];
            turboFFT_ZMUL(rPtr[14], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[30];
            turboFFT_ZMUL(rPtr[30], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[1];
            turboFFT_ZMUL(rPtr[1], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[17];
            turboFFT_ZMUL(rPtr[17], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[9];
            turboFFT_ZMUL(rPtr[9], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[25];
            turboFFT_ZMUL(rPtr[25], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[5];
            turboFFT_ZMUL(rPtr[5], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[21];
            turboFFT_ZMUL(rPtr[21], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[13];
            turboFFT_ZMUL(rPtr[13], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[29];
            turboFFT_ZMUL(rPtr[29], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[3];
            turboFFT_ZMUL(rPtr[3], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[19];
            turboFFT_ZMUL(rPtr[19], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[11];
            turboFFT_ZMUL(rPtr[11], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[27];
            turboFFT_ZMUL(rPtr[27], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[7];
            turboFFT_ZMUL(rPtr[7], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[23];
            turboFFT_ZMUL(rPtr[23], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[15];
            turboFFT_ZMUL(rPtr[15], tmp, angle);
            
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        
            tmp = rPtr[31];
            turboFFT_ZMUL(rPtr[31], tmp, angle);
            
            *(gPtr + 0) = rPtr[0];
            rPtr_4[0].x += rPtr[0].x;
            rPtr_4[0].y += rPtr[0].y;
            
            *(gPtr + 32768) = rPtr[16];
            rPtr_4[1].x += rPtr[16].x;
            rPtr_4[1].y += rPtr[16].y;
            
            *(gPtr + 65536) = rPtr[8];
            rPtr_4[2].x += rPtr[8].x;
            rPtr_4[2].y += rPtr[8].y;
            
            *(gPtr + 98304) = rPtr[24];
            rPtr_4[3].x += rPtr[24].x;
            rPtr_4[3].y += rPtr[24].y;
            
            *(gPtr + 131072) = rPtr[4];
            rPtr_4[4].x += rPtr[4].x;
            rPtr_4[4].y += rPtr[4].y;
            
            *(gPtr + 163840) = rPtr[20];
            rPtr_4[5].x += rPtr[20].x;
            rPtr_4[5].y += rPtr[20].y;
            
            *(gPtr + 196608) = rPtr[12];
            rPtr_4[6].x += rPtr[12].x;
            rPtr_4[6].y += rPtr[12].y;
            
            *(gPtr + 229376) = rPtr[28];
            rPtr_4[7].x += rPtr[28].x;
            rPtr_4[7].y += rPtr[28].y;
            
            *(gPtr + 262144) = rPtr[2];
            rPtr_4[8].x += rPtr[2].x;
            rPtr_4[8].y += rPtr[2].y;
            
            *(gPtr + 294912) = rPtr[18];
            rPtr_4[9].x += rPtr[18].x;
            rPtr_4[9].y += rPtr[18].y;
            
            *(gPtr + 327680) = rPtr[10];
            rPtr_4[10].x += rPtr[10].x;
            rPtr_4[10].y += rPtr[10].y;
            
            *(gPtr + 360448) = rPtr[26];
            rPtr_4[11].x += rPtr[26].x;
            rPtr_4[11].y += rPtr[26].y;
            
            *(gPtr + 393216) = rPtr[6];
            rPtr_4[12].x += rPtr[6].x;
            rPtr_4[12].y += rPtr[6].y;
            
            *(gPtr + 425984) = rPtr[22];
            rPtr_4[13].x += rPtr[22].x;
            rPtr_4[13].y += rPtr[22].y;
            
            *(gPtr + 458752) = rPtr[14];
            rPtr_4[14].x += rPtr[14].x;
            rPtr_4[14].y += rPtr[14].y;
            
            *(gPtr + 491520) = rPtr[30];
            rPtr_4[15].x += rPtr[30].x;
            rPtr_4[15].y += rPtr[30].y;
            
            *(gPtr + 524288) = rPtr[1];
            rPtr_4[16].x += rPtr[1].x;
            rPtr_4[16].y += rPtr[1].y;
            
            *(gPtr + 557056) = rPtr[17];
            rPtr_4[17].x += rPtr[17].x;
            rPtr_4[17].y += rPtr[17].y;
            
            *(gPtr + 589824) = rPtr[9];
            rPtr_4[18].x += rPtr[9].x;
            rPtr_4[18].y += rPtr[9].y;
            
            *(gPtr + 622592) = rPtr[25];
            rPtr_4[19].x += rPtr[25].x;
            rPtr_4[19].y += rPtr[25].y;
            
            *(gPtr + 655360) = rPtr[5];
            rPtr_4[20].x += rPtr[5].x;
            rPtr_4[20].y += rPtr[5].y;
            
            *(gPtr + 688128) = rPtr[21];
            rPtr_4[21].x += rPtr[21].x;
            rPtr_4[21].y += rPtr[21].y;
            
            *(gPtr + 720896) = rPtr[13];
            rPtr_4[22].x += rPtr[13].x;
            rPtr_4[22].y += rPtr[13].y;
            
            *(gPtr + 753664) = rPtr[29];
            rPtr_4[23].x += rPtr[29].x;
            rPtr_4[23].y += rPtr[29].y;
            
            *(gPtr + 786432) = rPtr[3];
            rPtr_4[24].x += rPtr[3].x;
            rPtr_4[24].y += rPtr[3].y;
            
            *(gPtr + 819200) = rPtr[19];
            rPtr_4[25].x += rPtr[19].x;
            rPtr_4[25].y += rPtr[19].y;
            
            *(gPtr + 851968) = rPtr[11];
            rPtr_4[26].x += rPtr[11].x;
            rPtr_4[26].y += rPtr[11].y;
            
            *(gPtr + 884736) = rPtr[27];
            rPtr_4[27].x += rPtr[27].x;
            rPtr_4[27].y += rPtr[27].y;
            
            *(gPtr + 917504) = rPtr[7];
            rPtr_4[28].x += rPtr[7].x;
            rPtr_4[28].y += rPtr[7].y;
            
            *(gPtr + 950272) = rPtr[23];
            rPtr_4[29].x += rPtr[23].x;
            rPtr_4[29].y += rPtr[23].y;
            
            *(gPtr + 983040) = rPtr[15];
            rPtr_4[30].x += rPtr[15].x;
            rPtr_4[30].y += rPtr[15].y;
            
            *(gPtr + 1015808) = rPtr[31];
            rPtr_4[31].x += rPtr[31].x;
            rPtr_4[31].y += rPtr[31].y;
            
    }
    
}
