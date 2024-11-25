
#include "../../../TurboFFT_radix_2_template.h"
template<>
__global__ void fft_radix_2<float2, 21, 1, 0, 0>(float2* inputs, float2* outputs, float2* twiddle, float2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    float2* shared = (float2*) ext_shared;
    int threadblock_per_SM = 4;
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
                bid_cnt < thread_bs && bid < (2097152 * BS + 8192 - 1) / 8192; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 4 * 1;
    
    gPtr += (bx % 1) * 2048 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 2048;
    bx = bx / 256;
    
    gPtr += tx % 4 * 2048;
    
    gPtr += (bx % BS * 2097152);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        rPtr[1] = *(gPtr + 64);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
        rPtr[2] = *(gPtr + 128);
        rPtr_3[2].x += rPtr[2].x;
        rPtr_3[2].y += rPtr[2].y;
        
        rPtr[3] = *(gPtr + 192);
        rPtr_3[3].x += rPtr[3].x;
        rPtr_3[3].y += rPtr[3].y;
        
        rPtr[4] = *(gPtr + 256);
        rPtr_3[4].x += rPtr[4].x;
        rPtr_3[4].y += rPtr[4].y;
        
        rPtr[5] = *(gPtr + 320);
        rPtr_3[5].x += rPtr[5].x;
        rPtr_3[5].y += rPtr[5].y;
        
        rPtr[6] = *(gPtr + 384);
        rPtr_3[6].x += rPtr[6].x;
        rPtr_3[6].y += rPtr[6].y;
        
        rPtr[7] = *(gPtr + 448);
        rPtr_3[7].x += rPtr[7].x;
        rPtr_3[7].y += rPtr[7].y;
        
        rPtr[8] = *(gPtr + 512);
        rPtr_3[8].x += rPtr[8].x;
        rPtr_3[8].y += rPtr[8].y;
        
        rPtr[9] = *(gPtr + 576);
        rPtr_3[9].x += rPtr[9].x;
        rPtr_3[9].y += rPtr[9].y;
        
        rPtr[10] = *(gPtr + 640);
        rPtr_3[10].x += rPtr[10].x;
        rPtr_3[10].y += rPtr[10].y;
        
        rPtr[11] = *(gPtr + 704);
        rPtr_3[11].x += rPtr[11].x;
        rPtr_3[11].y += rPtr[11].y;
        
        rPtr[12] = *(gPtr + 768);
        rPtr_3[12].x += rPtr[12].x;
        rPtr_3[12].y += rPtr[12].y;
        
        rPtr[13] = *(gPtr + 832);
        rPtr_3[13].x += rPtr[13].x;
        rPtr_3[13].y += rPtr[13].y;
        
        rPtr[14] = *(gPtr + 896);
        rPtr_3[14].x += rPtr[14].x;
        rPtr_3[14].y += rPtr[14].y;
        
        rPtr[15] = *(gPtr + 960);
        rPtr_3[15].x += rPtr[15].x;
        rPtr_3[15].y += rPtr[15].y;
        
        rPtr[16] = *(gPtr + 1024);
        rPtr_3[16].x += rPtr[16].x;
        rPtr_3[16].y += rPtr[16].y;
        
        rPtr[17] = *(gPtr + 1088);
        rPtr_3[17].x += rPtr[17].x;
        rPtr_3[17].y += rPtr[17].y;
        
        rPtr[18] = *(gPtr + 1152);
        rPtr_3[18].x += rPtr[18].x;
        rPtr_3[18].y += rPtr[18].y;
        
        rPtr[19] = *(gPtr + 1216);
        rPtr_3[19].x += rPtr[19].x;
        rPtr_3[19].y += rPtr[19].y;
        
        rPtr[20] = *(gPtr + 1280);
        rPtr_3[20].x += rPtr[20].x;
        rPtr_3[20].y += rPtr[20].y;
        
        rPtr[21] = *(gPtr + 1344);
        rPtr_3[21].x += rPtr[21].x;
        rPtr_3[21].y += rPtr[21].y;
        
        rPtr[22] = *(gPtr + 1408);
        rPtr_3[22].x += rPtr[22].x;
        rPtr_3[22].y += rPtr[22].y;
        
        rPtr[23] = *(gPtr + 1472);
        rPtr_3[23].x += rPtr[23].x;
        rPtr_3[23].y += rPtr[23].y;
        
        rPtr[24] = *(gPtr + 1536);
        rPtr_3[24].x += rPtr[24].x;
        rPtr_3[24].y += rPtr[24].y;
        
        rPtr[25] = *(gPtr + 1600);
        rPtr_3[25].x += rPtr[25].x;
        rPtr_3[25].y += rPtr[25].y;
        
        rPtr[26] = *(gPtr + 1664);
        rPtr_3[26].x += rPtr[26].x;
        rPtr_3[26].y += rPtr[26].y;
        
        rPtr[27] = *(gPtr + 1728);
        rPtr_3[27].x += rPtr[27].x;
        rPtr_3[27].y += rPtr[27].y;
        
        rPtr[28] = *(gPtr + 1792);
        rPtr_3[28].x += rPtr[28].x;
        rPtr_3[28].y += rPtr[28].y;
        
        rPtr[29] = *(gPtr + 1856);
        rPtr_3[29].x += rPtr[29].x;
        rPtr_3[29].y += rPtr[29].y;
        
        rPtr[30] = *(gPtr + 1920);
        rPtr_3[30].x += rPtr[30].x;
        rPtr_3[30].y += rPtr[30].y;
        
        rPtr[31] = *(gPtr + 1984);
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
    
    offset += ((tx / 1) % 4) * 1;
    
    j = tx / 4;
    
    offset += ((tx / 4) % 2) * 128;
    
    offset += ((tx / 8) % 32) * 256;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.003067961661145091f);
    delta_angle.y = __sinf(j * -0.003067961661145091f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 4] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 8] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 12] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 16] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 20] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 24] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 28] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 32] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 36] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 40] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 44] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 48] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 52] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 56] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 60] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 64] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 68] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 72] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 76] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 80] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 84] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 88] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 92] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 96] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 100] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 104] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 108] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 112] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 116] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 120] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 124] = rPtr[31];
    
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
    
    j = 0;
    offset  = 0;
    
    offset += ((tx / 1) % 4) * 1;
    
    offset += ((tx / 4) % 32) * 4;
    
    j = tx / 128;
    
    offset += ((tx / 128) % 2) * 4096;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.09817477315664291f);
    delta_angle.y = __sinf(j * -0.09817477315664291f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 128] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 256] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 384] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 512] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 640] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 768] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 896] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 1024] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 1152] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 1280] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 1408] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 1536] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 1664] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 1792] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 1920] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 2048] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 2176] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 2304] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 2432] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 2560] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 2688] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 2816] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 2944] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 3072] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 3200] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 3328] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 3456] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 3584] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 3712] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 3840] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 3968] = rPtr[31];
    
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
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
            
    bx = bid;
    tx = threadIdx.x;
    gPtr = outputs;
    
    gPtr += tx / 4 * 1024;
    
    gPtr += (bx % 1) * 2048 * 1024;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 1;
    bx = bx / 256;
    
    gPtr += tx % 4 * 1;
    
    gPtr += (bx % BS * 2097152);
    
            *(gPtr + 0) = rPtr[0];
            rPtr_4[0].x += rPtr[0].x;
            rPtr_4[0].y += rPtr[0].y;
            
            *(gPtr + 65536) = rPtr[1];
            rPtr_4[1].x += rPtr[1].x;
            rPtr_4[1].y += rPtr[1].y;
            
            *(gPtr + 131072) = rPtr[2];
            rPtr_4[2].x += rPtr[2].x;
            rPtr_4[2].y += rPtr[2].y;
            
            *(gPtr + 196608) = rPtr[3];
            rPtr_4[3].x += rPtr[3].x;
            rPtr_4[3].y += rPtr[3].y;
            
            *(gPtr + 262144) = rPtr[4];
            rPtr_4[4].x += rPtr[4].x;
            rPtr_4[4].y += rPtr[4].y;
            
            *(gPtr + 327680) = rPtr[5];
            rPtr_4[5].x += rPtr[5].x;
            rPtr_4[5].y += rPtr[5].y;
            
            *(gPtr + 393216) = rPtr[6];
            rPtr_4[6].x += rPtr[6].x;
            rPtr_4[6].y += rPtr[6].y;
            
            *(gPtr + 458752) = rPtr[7];
            rPtr_4[7].x += rPtr[7].x;
            rPtr_4[7].y += rPtr[7].y;
            
            *(gPtr + 524288) = rPtr[8];
            rPtr_4[8].x += rPtr[8].x;
            rPtr_4[8].y += rPtr[8].y;
            
            *(gPtr + 589824) = rPtr[9];
            rPtr_4[9].x += rPtr[9].x;
            rPtr_4[9].y += rPtr[9].y;
            
            *(gPtr + 655360) = rPtr[10];
            rPtr_4[10].x += rPtr[10].x;
            rPtr_4[10].y += rPtr[10].y;
            
            *(gPtr + 720896) = rPtr[11];
            rPtr_4[11].x += rPtr[11].x;
            rPtr_4[11].y += rPtr[11].y;
            
            *(gPtr + 786432) = rPtr[12];
            rPtr_4[12].x += rPtr[12].x;
            rPtr_4[12].y += rPtr[12].y;
            
            *(gPtr + 851968) = rPtr[13];
            rPtr_4[13].x += rPtr[13].x;
            rPtr_4[13].y += rPtr[13].y;
            
            *(gPtr + 917504) = rPtr[14];
            rPtr_4[14].x += rPtr[14].x;
            rPtr_4[14].y += rPtr[14].y;
            
            *(gPtr + 983040) = rPtr[15];
            rPtr_4[15].x += rPtr[15].x;
            rPtr_4[15].y += rPtr[15].y;
            
            *(gPtr + 1048576) = rPtr[16];
            rPtr_4[16].x += rPtr[16].x;
            rPtr_4[16].y += rPtr[16].y;
            
            *(gPtr + 1114112) = rPtr[17];
            rPtr_4[17].x += rPtr[17].x;
            rPtr_4[17].y += rPtr[17].y;
            
            *(gPtr + 1179648) = rPtr[18];
            rPtr_4[18].x += rPtr[18].x;
            rPtr_4[18].y += rPtr[18].y;
            
            *(gPtr + 1245184) = rPtr[19];
            rPtr_4[19].x += rPtr[19].x;
            rPtr_4[19].y += rPtr[19].y;
            
            *(gPtr + 1310720) = rPtr[20];
            rPtr_4[20].x += rPtr[20].x;
            rPtr_4[20].y += rPtr[20].y;
            
            *(gPtr + 1376256) = rPtr[21];
            rPtr_4[21].x += rPtr[21].x;
            rPtr_4[21].y += rPtr[21].y;
            
            *(gPtr + 1441792) = rPtr[22];
            rPtr_4[22].x += rPtr[22].x;
            rPtr_4[22].y += rPtr[22].y;
            
            *(gPtr + 1507328) = rPtr[23];
            rPtr_4[23].x += rPtr[23].x;
            rPtr_4[23].y += rPtr[23].y;
            
            *(gPtr + 1572864) = rPtr[24];
            rPtr_4[24].x += rPtr[24].x;
            rPtr_4[24].y += rPtr[24].y;
            
            *(gPtr + 1638400) = rPtr[25];
            rPtr_4[25].x += rPtr[25].x;
            rPtr_4[25].y += rPtr[25].y;
            
            *(gPtr + 1703936) = rPtr[26];
            rPtr_4[26].x += rPtr[26].x;
            rPtr_4[26].y += rPtr[26].y;
            
            *(gPtr + 1769472) = rPtr[27];
            rPtr_4[27].x += rPtr[27].x;
            rPtr_4[27].y += rPtr[27].y;
            
            *(gPtr + 1835008) = rPtr[28];
            rPtr_4[28].x += rPtr[28].x;
            rPtr_4[28].y += rPtr[28].y;
            
            *(gPtr + 1900544) = rPtr[29];
            rPtr_4[29].x += rPtr[29].x;
            rPtr_4[29].y += rPtr[29].y;
            
            *(gPtr + 1966080) = rPtr[30];
            rPtr_4[30].x += rPtr[30].x;
            rPtr_4[30].y += rPtr[30].y;
            
            *(gPtr + 2031616) = rPtr[31];
            rPtr_4[31].x += rPtr[31].x;
            rPtr_4[31].y += rPtr[31].y;
            
    }
    
}

#include "../../../TurboFFT_radix_2_template.h"
template<>
__global__ void fft_radix_2<float2, 21, 1, 1, 0>(float2* inputs, float2* outputs, float2* twiddle, float2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    float2* shared = (float2*) ext_shared;
    int threadblock_per_SM = 4;
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
    
    rPtr_2[0] = *(checksum_DFT + 2048 - 2 + tx + 0);
    shPtr[tx + 0] = rPtr_2[0];
    
    rPtr_2[1] = *(checksum_DFT + 2048 - 2 + tx + 256);
    shPtr[tx + 256] = rPtr_2[1];
    
    rPtr_2[2] = *(checksum_DFT + 2048 - 2 + tx + 512);
    shPtr[tx + 512] = rPtr_2[2];
    
    rPtr_2[3] = *(checksum_DFT + 2048 - 2 + tx + 768);
    shPtr[tx + 768] = rPtr_2[3];
    
    rPtr_2[4] = *(checksum_DFT + 2048 - 2 + tx + 1024);
    shPtr[tx + 1024] = rPtr_2[4];
    
    rPtr_2[5] = *(checksum_DFT + 2048 - 2 + tx + 1280);
    shPtr[tx + 1280] = rPtr_2[5];
    
    rPtr_2[6] = *(checksum_DFT + 2048 - 2 + tx + 1536);
    shPtr[tx + 1536] = rPtr_2[6];
    
    rPtr_2[7] = *(checksum_DFT + 2048 - 2 + tx + 1792);
    shPtr[tx + 1792] = rPtr_2[7];
    
    __syncthreads();
    tmp_1.x = 0;
    tmp_1.y = 0;
    tmp_2.x = 0;
    tmp_2.y = 0;
    tmp_3.x = 0;
    tmp_3.y = 0;
    
    rPtr_2[0] = *(shPtr +  tx / 4 + 0);
    rPtr_3[0].x = 0; rPtr_3[0].y = 0;
    rPtr_4[0].x = 0; rPtr_4[0].y = 0;
    
    rPtr_2[1] = *(shPtr +  tx / 4 + 64);
    rPtr_3[1].x = 0; rPtr_3[1].y = 0;
    rPtr_4[1].x = 0; rPtr_4[1].y = 0;
    
    rPtr_2[2] = *(shPtr +  tx / 4 + 128);
    rPtr_3[2].x = 0; rPtr_3[2].y = 0;
    rPtr_4[2].x = 0; rPtr_4[2].y = 0;
    
    rPtr_2[3] = *(shPtr +  tx / 4 + 192);
    rPtr_3[3].x = 0; rPtr_3[3].y = 0;
    rPtr_4[3].x = 0; rPtr_4[3].y = 0;
    
    rPtr_2[4] = *(shPtr +  tx / 4 + 256);
    rPtr_3[4].x = 0; rPtr_3[4].y = 0;
    rPtr_4[4].x = 0; rPtr_4[4].y = 0;
    
    rPtr_2[5] = *(shPtr +  tx / 4 + 320);
    rPtr_3[5].x = 0; rPtr_3[5].y = 0;
    rPtr_4[5].x = 0; rPtr_4[5].y = 0;
    
    rPtr_2[6] = *(shPtr +  tx / 4 + 384);
    rPtr_3[6].x = 0; rPtr_3[6].y = 0;
    rPtr_4[6].x = 0; rPtr_4[6].y = 0;
    
    rPtr_2[7] = *(shPtr +  tx / 4 + 448);
    rPtr_3[7].x = 0; rPtr_3[7].y = 0;
    rPtr_4[7].x = 0; rPtr_4[7].y = 0;
    
    rPtr_2[8] = *(shPtr +  tx / 4 + 512);
    rPtr_3[8].x = 0; rPtr_3[8].y = 0;
    rPtr_4[8].x = 0; rPtr_4[8].y = 0;
    
    rPtr_2[9] = *(shPtr +  tx / 4 + 576);
    rPtr_3[9].x = 0; rPtr_3[9].y = 0;
    rPtr_4[9].x = 0; rPtr_4[9].y = 0;
    
    rPtr_2[10] = *(shPtr +  tx / 4 + 640);
    rPtr_3[10].x = 0; rPtr_3[10].y = 0;
    rPtr_4[10].x = 0; rPtr_4[10].y = 0;
    
    rPtr_2[11] = *(shPtr +  tx / 4 + 704);
    rPtr_3[11].x = 0; rPtr_3[11].y = 0;
    rPtr_4[11].x = 0; rPtr_4[11].y = 0;
    
    rPtr_2[12] = *(shPtr +  tx / 4 + 768);
    rPtr_3[12].x = 0; rPtr_3[12].y = 0;
    rPtr_4[12].x = 0; rPtr_4[12].y = 0;
    
    rPtr_2[13] = *(shPtr +  tx / 4 + 832);
    rPtr_3[13].x = 0; rPtr_3[13].y = 0;
    rPtr_4[13].x = 0; rPtr_4[13].y = 0;
    
    rPtr_2[14] = *(shPtr +  tx / 4 + 896);
    rPtr_3[14].x = 0; rPtr_3[14].y = 0;
    rPtr_4[14].x = 0; rPtr_4[14].y = 0;
    
    rPtr_2[15] = *(shPtr +  tx / 4 + 960);
    rPtr_3[15].x = 0; rPtr_3[15].y = 0;
    rPtr_4[15].x = 0; rPtr_4[15].y = 0;
    
    rPtr_2[16] = *(shPtr +  tx / 4 + 1024);
    rPtr_3[16].x = 0; rPtr_3[16].y = 0;
    rPtr_4[16].x = 0; rPtr_4[16].y = 0;
    
    rPtr_2[17] = *(shPtr +  tx / 4 + 1088);
    rPtr_3[17].x = 0; rPtr_3[17].y = 0;
    rPtr_4[17].x = 0; rPtr_4[17].y = 0;
    
    rPtr_2[18] = *(shPtr +  tx / 4 + 1152);
    rPtr_3[18].x = 0; rPtr_3[18].y = 0;
    rPtr_4[18].x = 0; rPtr_4[18].y = 0;
    
    rPtr_2[19] = *(shPtr +  tx / 4 + 1216);
    rPtr_3[19].x = 0; rPtr_3[19].y = 0;
    rPtr_4[19].x = 0; rPtr_4[19].y = 0;
    
    rPtr_2[20] = *(shPtr +  tx / 4 + 1280);
    rPtr_3[20].x = 0; rPtr_3[20].y = 0;
    rPtr_4[20].x = 0; rPtr_4[20].y = 0;
    
    rPtr_2[21] = *(shPtr +  tx / 4 + 1344);
    rPtr_3[21].x = 0; rPtr_3[21].y = 0;
    rPtr_4[21].x = 0; rPtr_4[21].y = 0;
    
    rPtr_2[22] = *(shPtr +  tx / 4 + 1408);
    rPtr_3[22].x = 0; rPtr_3[22].y = 0;
    rPtr_4[22].x = 0; rPtr_4[22].y = 0;
    
    rPtr_2[23] = *(shPtr +  tx / 4 + 1472);
    rPtr_3[23].x = 0; rPtr_3[23].y = 0;
    rPtr_4[23].x = 0; rPtr_4[23].y = 0;
    
    rPtr_2[24] = *(shPtr +  tx / 4 + 1536);
    rPtr_3[24].x = 0; rPtr_3[24].y = 0;
    rPtr_4[24].x = 0; rPtr_4[24].y = 0;
    
    rPtr_2[25] = *(shPtr +  tx / 4 + 1600);
    rPtr_3[25].x = 0; rPtr_3[25].y = 0;
    rPtr_4[25].x = 0; rPtr_4[25].y = 0;
    
    rPtr_2[26] = *(shPtr +  tx / 4 + 1664);
    rPtr_3[26].x = 0; rPtr_3[26].y = 0;
    rPtr_4[26].x = 0; rPtr_4[26].y = 0;
    
    rPtr_2[27] = *(shPtr +  tx / 4 + 1728);
    rPtr_3[27].x = 0; rPtr_3[27].y = 0;
    rPtr_4[27].x = 0; rPtr_4[27].y = 0;
    
    rPtr_2[28] = *(shPtr +  tx / 4 + 1792);
    rPtr_3[28].x = 0; rPtr_3[28].y = 0;
    rPtr_4[28].x = 0; rPtr_4[28].y = 0;
    
    rPtr_2[29] = *(shPtr +  tx / 4 + 1856);
    rPtr_3[29].x = 0; rPtr_3[29].y = 0;
    rPtr_4[29].x = 0; rPtr_4[29].y = 0;
    
    rPtr_2[30] = *(shPtr +  tx / 4 + 1920);
    rPtr_3[30].x = 0; rPtr_3[30].y = 0;
    rPtr_4[30].x = 0; rPtr_4[30].y = 0;
    
    rPtr_2[31] = *(shPtr +  tx / 4 + 1984);
    rPtr_3[31].x = 0; rPtr_3[31].y = 0;
    rPtr_4[31].x = 0; rPtr_4[31].y = 0;
    
    __syncthreads();
    int bid = 0;
    for(bid = (blockIdx.x / tb_gap) * tb_gap * thread_bs + blockIdx.x % tb_gap;
                bid_cnt < thread_bs && bid < (2097152 * BS + 8192 - 1) / 8192; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 4 * 1;
    
    gPtr += (bx % 1) * 2048 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 2048;
    bx = bx / 256;
    
    gPtr += tx % 4 * 2048;
    
    gPtr += (bx % BS * 2097152);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        // tmp = checksum_DFT[tx / 4 + 0];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[0], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[0], rPtr_2[0])
        turboFFT_ZMUL(tmp, rPtr[0], rPtr_2[0])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[1] = *(gPtr + 64);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
        // tmp = checksum_DFT[tx / 4 + 64];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[1], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[1], rPtr_2[1])
        turboFFT_ZMUL(tmp, rPtr[1], rPtr_2[1])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[2] = *(gPtr + 128);
        rPtr_3[2].x += rPtr[2].x;
        rPtr_3[2].y += rPtr[2].y;
        
        // tmp = checksum_DFT[tx / 4 + 128];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[2], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[2], rPtr_2[2])
        turboFFT_ZMUL(tmp, rPtr[2], rPtr_2[2])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[3] = *(gPtr + 192);
        rPtr_3[3].x += rPtr[3].x;
        rPtr_3[3].y += rPtr[3].y;
        
        // tmp = checksum_DFT[tx / 4 + 192];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[3], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[3], rPtr_2[3])
        turboFFT_ZMUL(tmp, rPtr[3], rPtr_2[3])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[4] = *(gPtr + 256);
        rPtr_3[4].x += rPtr[4].x;
        rPtr_3[4].y += rPtr[4].y;
        
        // tmp = checksum_DFT[tx / 4 + 256];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[4], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[4], rPtr_2[4])
        turboFFT_ZMUL(tmp, rPtr[4], rPtr_2[4])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[5] = *(gPtr + 320);
        rPtr_3[5].x += rPtr[5].x;
        rPtr_3[5].y += rPtr[5].y;
        
        // tmp = checksum_DFT[tx / 4 + 320];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[5], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[5], rPtr_2[5])
        turboFFT_ZMUL(tmp, rPtr[5], rPtr_2[5])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[6] = *(gPtr + 384);
        rPtr_3[6].x += rPtr[6].x;
        rPtr_3[6].y += rPtr[6].y;
        
        // tmp = checksum_DFT[tx / 4 + 384];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[6], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[6], rPtr_2[6])
        turboFFT_ZMUL(tmp, rPtr[6], rPtr_2[6])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[7] = *(gPtr + 448);
        rPtr_3[7].x += rPtr[7].x;
        rPtr_3[7].y += rPtr[7].y;
        
        // tmp = checksum_DFT[tx / 4 + 448];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[7], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[7], rPtr_2[7])
        turboFFT_ZMUL(tmp, rPtr[7], rPtr_2[7])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[8] = *(gPtr + 512);
        rPtr_3[8].x += rPtr[8].x;
        rPtr_3[8].y += rPtr[8].y;
        
        // tmp = checksum_DFT[tx / 4 + 512];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[8], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[8], rPtr_2[8])
        turboFFT_ZMUL(tmp, rPtr[8], rPtr_2[8])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[9] = *(gPtr + 576);
        rPtr_3[9].x += rPtr[9].x;
        rPtr_3[9].y += rPtr[9].y;
        
        // tmp = checksum_DFT[tx / 4 + 576];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[9], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[9], rPtr_2[9])
        turboFFT_ZMUL(tmp, rPtr[9], rPtr_2[9])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[10] = *(gPtr + 640);
        rPtr_3[10].x += rPtr[10].x;
        rPtr_3[10].y += rPtr[10].y;
        
        // tmp = checksum_DFT[tx / 4 + 640];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[10], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[10], rPtr_2[10])
        turboFFT_ZMUL(tmp, rPtr[10], rPtr_2[10])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[11] = *(gPtr + 704);
        rPtr_3[11].x += rPtr[11].x;
        rPtr_3[11].y += rPtr[11].y;
        
        // tmp = checksum_DFT[tx / 4 + 704];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[11], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[11], rPtr_2[11])
        turboFFT_ZMUL(tmp, rPtr[11], rPtr_2[11])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[12] = *(gPtr + 768);
        rPtr_3[12].x += rPtr[12].x;
        rPtr_3[12].y += rPtr[12].y;
        
        // tmp = checksum_DFT[tx / 4 + 768];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[12], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[12], rPtr_2[12])
        turboFFT_ZMUL(tmp, rPtr[12], rPtr_2[12])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[13] = *(gPtr + 832);
        rPtr_3[13].x += rPtr[13].x;
        rPtr_3[13].y += rPtr[13].y;
        
        // tmp = checksum_DFT[tx / 4 + 832];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[13], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[13], rPtr_2[13])
        turboFFT_ZMUL(tmp, rPtr[13], rPtr_2[13])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[14] = *(gPtr + 896);
        rPtr_3[14].x += rPtr[14].x;
        rPtr_3[14].y += rPtr[14].y;
        
        // tmp = checksum_DFT[tx / 4 + 896];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[14], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[14], rPtr_2[14])
        turboFFT_ZMUL(tmp, rPtr[14], rPtr_2[14])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[15] = *(gPtr + 960);
        rPtr_3[15].x += rPtr[15].x;
        rPtr_3[15].y += rPtr[15].y;
        
        // tmp = checksum_DFT[tx / 4 + 960];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[15], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[15], rPtr_2[15])
        turboFFT_ZMUL(tmp, rPtr[15], rPtr_2[15])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[16] = *(gPtr + 1024);
        rPtr_3[16].x += rPtr[16].x;
        rPtr_3[16].y += rPtr[16].y;
        
        // tmp = checksum_DFT[tx / 4 + 1024];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[16], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[16], rPtr_2[16])
        turboFFT_ZMUL(tmp, rPtr[16], rPtr_2[16])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[17] = *(gPtr + 1088);
        rPtr_3[17].x += rPtr[17].x;
        rPtr_3[17].y += rPtr[17].y;
        
        // tmp = checksum_DFT[tx / 4 + 1088];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[17], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[17], rPtr_2[17])
        turboFFT_ZMUL(tmp, rPtr[17], rPtr_2[17])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[18] = *(gPtr + 1152);
        rPtr_3[18].x += rPtr[18].x;
        rPtr_3[18].y += rPtr[18].y;
        
        // tmp = checksum_DFT[tx / 4 + 1152];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[18], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[18], rPtr_2[18])
        turboFFT_ZMUL(tmp, rPtr[18], rPtr_2[18])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[19] = *(gPtr + 1216);
        rPtr_3[19].x += rPtr[19].x;
        rPtr_3[19].y += rPtr[19].y;
        
        // tmp = checksum_DFT[tx / 4 + 1216];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[19], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[19], rPtr_2[19])
        turboFFT_ZMUL(tmp, rPtr[19], rPtr_2[19])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[20] = *(gPtr + 1280);
        rPtr_3[20].x += rPtr[20].x;
        rPtr_3[20].y += rPtr[20].y;
        
        // tmp = checksum_DFT[tx / 4 + 1280];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[20], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[20], rPtr_2[20])
        turboFFT_ZMUL(tmp, rPtr[20], rPtr_2[20])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[21] = *(gPtr + 1344);
        rPtr_3[21].x += rPtr[21].x;
        rPtr_3[21].y += rPtr[21].y;
        
        // tmp = checksum_DFT[tx / 4 + 1344];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[21], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[21], rPtr_2[21])
        turboFFT_ZMUL(tmp, rPtr[21], rPtr_2[21])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[22] = *(gPtr + 1408);
        rPtr_3[22].x += rPtr[22].x;
        rPtr_3[22].y += rPtr[22].y;
        
        // tmp = checksum_DFT[tx / 4 + 1408];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[22], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[22], rPtr_2[22])
        turboFFT_ZMUL(tmp, rPtr[22], rPtr_2[22])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[23] = *(gPtr + 1472);
        rPtr_3[23].x += rPtr[23].x;
        rPtr_3[23].y += rPtr[23].y;
        
        // tmp = checksum_DFT[tx / 4 + 1472];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[23], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[23], rPtr_2[23])
        turboFFT_ZMUL(tmp, rPtr[23], rPtr_2[23])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[24] = *(gPtr + 1536);
        rPtr_3[24].x += rPtr[24].x;
        rPtr_3[24].y += rPtr[24].y;
        
        // tmp = checksum_DFT[tx / 4 + 1536];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[24], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[24], rPtr_2[24])
        turboFFT_ZMUL(tmp, rPtr[24], rPtr_2[24])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[25] = *(gPtr + 1600);
        rPtr_3[25].x += rPtr[25].x;
        rPtr_3[25].y += rPtr[25].y;
        
        // tmp = checksum_DFT[tx / 4 + 1600];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[25], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[25], rPtr_2[25])
        turboFFT_ZMUL(tmp, rPtr[25], rPtr_2[25])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[26] = *(gPtr + 1664);
        rPtr_3[26].x += rPtr[26].x;
        rPtr_3[26].y += rPtr[26].y;
        
        // tmp = checksum_DFT[tx / 4 + 1664];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[26], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[26], rPtr_2[26])
        turboFFT_ZMUL(tmp, rPtr[26], rPtr_2[26])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[27] = *(gPtr + 1728);
        rPtr_3[27].x += rPtr[27].x;
        rPtr_3[27].y += rPtr[27].y;
        
        // tmp = checksum_DFT[tx / 4 + 1728];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[27], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[27], rPtr_2[27])
        turboFFT_ZMUL(tmp, rPtr[27], rPtr_2[27])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[28] = *(gPtr + 1792);
        rPtr_3[28].x += rPtr[28].x;
        rPtr_3[28].y += rPtr[28].y;
        
        // tmp = checksum_DFT[tx / 4 + 1792];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[28], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[28], rPtr_2[28])
        turboFFT_ZMUL(tmp, rPtr[28], rPtr_2[28])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[29] = *(gPtr + 1856);
        rPtr_3[29].x += rPtr[29].x;
        rPtr_3[29].y += rPtr[29].y;
        
        // tmp = checksum_DFT[tx / 4 + 1856];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[29], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[29], rPtr_2[29])
        turboFFT_ZMUL(tmp, rPtr[29], rPtr_2[29])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[30] = *(gPtr + 1920);
        rPtr_3[30].x += rPtr[30].x;
        rPtr_3[30].y += rPtr[30].y;
        
        // tmp = checksum_DFT[tx / 4 + 1920];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[30], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[30], rPtr_2[30])
        turboFFT_ZMUL(tmp, rPtr[30], rPtr_2[30])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[31] = *(gPtr + 1984);
        rPtr_3[31].x += rPtr[31].x;
        rPtr_3[31].y += rPtr[31].y;
        
        // tmp = checksum_DFT[tx / 4 + 1984];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[31], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[31], rPtr_2[31])
        turboFFT_ZMUL(tmp, rPtr[31], rPtr_2[31])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        // tmp_3.x += bid_cnt * (rPtr[0].x + rPtr[0].y) * 2048;
        
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
    
    offset += ((tx / 1) % 4) * 1;
    
    j = tx / 4;
    
    offset += ((tx / 4) % 2) * 128;
    
    offset += ((tx / 8) % 32) * 256;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.003067961661145091f);
    delta_angle.y = __sinf(j * -0.003067961661145091f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 4] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 8] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 12] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 16] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 20] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 24] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 28] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 32] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 36] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 40] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 44] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 48] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 52] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 56] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 60] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 64] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 68] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 72] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 76] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 80] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 84] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 88] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 92] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 96] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 100] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 104] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 108] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 112] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 116] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 120] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 124] = rPtr[31];
    
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
    
    j = 0;
    offset  = 0;
    
    offset += ((tx / 1) % 4) * 1;
    
    offset += ((tx / 4) % 32) * 4;
    
    j = tx / 128;
    
    offset += ((tx / 128) % 2) * 4096;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.09817477315664291f);
    delta_angle.y = __sinf(j * -0.09817477315664291f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 128] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 256] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 384] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 512] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 640] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 768] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 896] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 1024] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 1152] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 1280] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 1408] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 1536] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 1664] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 1792] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 1920] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 2048] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 2176] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 2304] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 2432] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 2560] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 2688] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 2816] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 2944] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 3072] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 3200] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 3328] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 3456] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 3584] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 3712] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 3840] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 3968] = rPtr[31];
    
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
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
            
    bx = bid;
    tx = threadIdx.x;
    gPtr = outputs;
    
    gPtr += tx / 4 * 1024;
    
    gPtr += (bx % 1) * 2048 * 1024;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 1;
    bx = bx / 256;
    
    gPtr += tx % 4 * 1;
    
    gPtr += (bx % BS * 2097152);
    
        // 1's vector
        // tmp_3.y -=  (rPtr[0].y + rPtr[0].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[0].y + rPtr[0].x);
        turboFFT_ZMUL(tmp, rPtr[0],r[(0 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[0], r[(0 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[0], r[(0 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[1].y + rPtr[1].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[1].y + rPtr[1].x);
        turboFFT_ZMUL(tmp, rPtr[1],r[(64 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[1], r[(64 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[1], r[(64 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[2].y + rPtr[2].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[2].y + rPtr[2].x);
        turboFFT_ZMUL(tmp, rPtr[2],r[(128 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[2], r[(128 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[2], r[(128 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[3].y + rPtr[3].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[3].y + rPtr[3].x);
        turboFFT_ZMUL(tmp, rPtr[3],r[(192 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[3], r[(192 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[3], r[(192 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[4].y + rPtr[4].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[4].y + rPtr[4].x);
        turboFFT_ZMUL(tmp, rPtr[4],r[(256 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[4], r[(256 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[4], r[(256 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[5].y + rPtr[5].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[5].y + rPtr[5].x);
        turboFFT_ZMUL(tmp, rPtr[5],r[(320 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[5], r[(320 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[5], r[(320 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[6].y + rPtr[6].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[6].y + rPtr[6].x);
        turboFFT_ZMUL(tmp, rPtr[6],r[(384 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[6], r[(384 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[6], r[(384 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[7].y + rPtr[7].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[7].y + rPtr[7].x);
        turboFFT_ZMUL(tmp, rPtr[7],r[(448 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[7], r[(448 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[7], r[(448 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[8].y + rPtr[8].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[8].y + rPtr[8].x);
        turboFFT_ZMUL(tmp, rPtr[8],r[(512 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[8], r[(512 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[8], r[(512 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[9].y + rPtr[9].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[9].y + rPtr[9].x);
        turboFFT_ZMUL(tmp, rPtr[9],r[(576 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[9], r[(576 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[9], r[(576 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[10].y + rPtr[10].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[10].y + rPtr[10].x);
        turboFFT_ZMUL(tmp, rPtr[10],r[(640 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[10], r[(640 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[10], r[(640 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[11].y + rPtr[11].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[11].y + rPtr[11].x);
        turboFFT_ZMUL(tmp, rPtr[11],r[(704 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[11], r[(704 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[11], r[(704 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[12].y + rPtr[12].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[12].y + rPtr[12].x);
        turboFFT_ZMUL(tmp, rPtr[12],r[(768 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[12], r[(768 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[12], r[(768 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[13].y + rPtr[13].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[13].y + rPtr[13].x);
        turboFFT_ZMUL(tmp, rPtr[13],r[(832 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[13], r[(832 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[13], r[(832 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[14].y + rPtr[14].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[14].y + rPtr[14].x);
        turboFFT_ZMUL(tmp, rPtr[14],r[(896 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[14], r[(896 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[14], r[(896 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[15].y + rPtr[15].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[15].y + rPtr[15].x);
        turboFFT_ZMUL(tmp, rPtr[15],r[(960 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[15], r[(960 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[15], r[(960 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[16].y + rPtr[16].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[16].y + rPtr[16].x);
        turboFFT_ZMUL(tmp, rPtr[16],r[(1024 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[16], r[(1024 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[16], r[(1024 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[17].y + rPtr[17].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[17].y + rPtr[17].x);
        turboFFT_ZMUL(tmp, rPtr[17],r[(1088 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[17], r[(1088 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[17], r[(1088 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[18].y + rPtr[18].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[18].y + rPtr[18].x);
        turboFFT_ZMUL(tmp, rPtr[18],r[(1152 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[18], r[(1152 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[18], r[(1152 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[19].y + rPtr[19].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[19].y + rPtr[19].x);
        turboFFT_ZMUL(tmp, rPtr[19],r[(1216 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[19], r[(1216 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[19], r[(1216 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[20].y + rPtr[20].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[20].y + rPtr[20].x);
        turboFFT_ZMUL(tmp, rPtr[20],r[(1280 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[20], r[(1280 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[20], r[(1280 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[21].y + rPtr[21].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[21].y + rPtr[21].x);
        turboFFT_ZMUL(tmp, rPtr[21],r[(1344 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[21], r[(1344 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[21], r[(1344 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[22].y + rPtr[22].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[22].y + rPtr[22].x);
        turboFFT_ZMUL(tmp, rPtr[22],r[(1408 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[22], r[(1408 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[22], r[(1408 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[23].y + rPtr[23].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[23].y + rPtr[23].x);
        turboFFT_ZMUL(tmp, rPtr[23],r[(1472 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[23], r[(1472 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[23], r[(1472 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[24].y + rPtr[24].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[24].y + rPtr[24].x);
        turboFFT_ZMUL(tmp, rPtr[24],r[(1536 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[24], r[(1536 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[24], r[(1536 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[25].y + rPtr[25].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[25].y + rPtr[25].x);
        turboFFT_ZMUL(tmp, rPtr[25],r[(1600 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[25], r[(1600 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[25], r[(1600 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[26].y + rPtr[26].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[26].y + rPtr[26].x);
        turboFFT_ZMUL(tmp, rPtr[26],r[(1664 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[26], r[(1664 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[26], r[(1664 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[27].y + rPtr[27].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[27].y + rPtr[27].x);
        turboFFT_ZMUL(tmp, rPtr[27],r[(1728 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[27], r[(1728 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[27], r[(1728 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[28].y + rPtr[28].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[28].y + rPtr[28].x);
        turboFFT_ZMUL(tmp, rPtr[28],r[(1792 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[28], r[(1792 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[28], r[(1792 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[29].y + rPtr[29].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[29].y + rPtr[29].x);
        turboFFT_ZMUL(tmp, rPtr[29],r[(1856 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[29], r[(1856 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[29], r[(1856 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[30].y + rPtr[30].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[30].y + rPtr[30].x);
        turboFFT_ZMUL(tmp, rPtr[30],r[(1920 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[30], r[(1920 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[30], r[(1920 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[31].y + rPtr[31].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[31].y + rPtr[31].x);
        turboFFT_ZMUL(tmp, rPtr[31],r[(1984 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[31], r[(1984 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[31], r[(1984 + tx / 4) % 3])
        
            *(gPtr + 0) = rPtr[0];
            rPtr_4[0].x += rPtr[0].x;
            rPtr_4[0].y += rPtr[0].y;
            
            *(gPtr + 65536) = rPtr[1];
            rPtr_4[1].x += rPtr[1].x;
            rPtr_4[1].y += rPtr[1].y;
            
            *(gPtr + 131072) = rPtr[2];
            rPtr_4[2].x += rPtr[2].x;
            rPtr_4[2].y += rPtr[2].y;
            
            *(gPtr + 196608) = rPtr[3];
            rPtr_4[3].x += rPtr[3].x;
            rPtr_4[3].y += rPtr[3].y;
            
            *(gPtr + 262144) = rPtr[4];
            rPtr_4[4].x += rPtr[4].x;
            rPtr_4[4].y += rPtr[4].y;
            
            *(gPtr + 327680) = rPtr[5];
            rPtr_4[5].x += rPtr[5].x;
            rPtr_4[5].y += rPtr[5].y;
            
            *(gPtr + 393216) = rPtr[6];
            rPtr_4[6].x += rPtr[6].x;
            rPtr_4[6].y += rPtr[6].y;
            
            *(gPtr + 458752) = rPtr[7];
            rPtr_4[7].x += rPtr[7].x;
            rPtr_4[7].y += rPtr[7].y;
            
            *(gPtr + 524288) = rPtr[8];
            rPtr_4[8].x += rPtr[8].x;
            rPtr_4[8].y += rPtr[8].y;
            
            *(gPtr + 589824) = rPtr[9];
            rPtr_4[9].x += rPtr[9].x;
            rPtr_4[9].y += rPtr[9].y;
            
            *(gPtr + 655360) = rPtr[10];
            rPtr_4[10].x += rPtr[10].x;
            rPtr_4[10].y += rPtr[10].y;
            
            *(gPtr + 720896) = rPtr[11];
            rPtr_4[11].x += rPtr[11].x;
            rPtr_4[11].y += rPtr[11].y;
            
            *(gPtr + 786432) = rPtr[12];
            rPtr_4[12].x += rPtr[12].x;
            rPtr_4[12].y += rPtr[12].y;
            
            *(gPtr + 851968) = rPtr[13];
            rPtr_4[13].x += rPtr[13].x;
            rPtr_4[13].y += rPtr[13].y;
            
            *(gPtr + 917504) = rPtr[14];
            rPtr_4[14].x += rPtr[14].x;
            rPtr_4[14].y += rPtr[14].y;
            
            *(gPtr + 983040) = rPtr[15];
            rPtr_4[15].x += rPtr[15].x;
            rPtr_4[15].y += rPtr[15].y;
            
            *(gPtr + 1048576) = rPtr[16];
            rPtr_4[16].x += rPtr[16].x;
            rPtr_4[16].y += rPtr[16].y;
            
            *(gPtr + 1114112) = rPtr[17];
            rPtr_4[17].x += rPtr[17].x;
            rPtr_4[17].y += rPtr[17].y;
            
            *(gPtr + 1179648) = rPtr[18];
            rPtr_4[18].x += rPtr[18].x;
            rPtr_4[18].y += rPtr[18].y;
            
            *(gPtr + 1245184) = rPtr[19];
            rPtr_4[19].x += rPtr[19].x;
            rPtr_4[19].y += rPtr[19].y;
            
            *(gPtr + 1310720) = rPtr[20];
            rPtr_4[20].x += rPtr[20].x;
            rPtr_4[20].y += rPtr[20].y;
            
            *(gPtr + 1376256) = rPtr[21];
            rPtr_4[21].x += rPtr[21].x;
            rPtr_4[21].y += rPtr[21].y;
            
            *(gPtr + 1441792) = rPtr[22];
            rPtr_4[22].x += rPtr[22].x;
            rPtr_4[22].y += rPtr[22].y;
            
            *(gPtr + 1507328) = rPtr[23];
            rPtr_4[23].x += rPtr[23].x;
            rPtr_4[23].y += rPtr[23].y;
            
            *(gPtr + 1572864) = rPtr[24];
            rPtr_4[24].x += rPtr[24].x;
            rPtr_4[24].y += rPtr[24].y;
            
            *(gPtr + 1638400) = rPtr[25];
            rPtr_4[25].x += rPtr[25].x;
            rPtr_4[25].y += rPtr[25].y;
            
            *(gPtr + 1703936) = rPtr[26];
            rPtr_4[26].x += rPtr[26].x;
            rPtr_4[26].y += rPtr[26].y;
            
            *(gPtr + 1769472) = rPtr[27];
            rPtr_4[27].x += rPtr[27].x;
            rPtr_4[27].y += rPtr[27].y;
            
            *(gPtr + 1835008) = rPtr[28];
            rPtr_4[28].x += rPtr[28].x;
            rPtr_4[28].y += rPtr[28].y;
            
            *(gPtr + 1900544) = rPtr[29];
            rPtr_4[29].x += rPtr[29].x;
            rPtr_4[29].y += rPtr[29].y;
            
            *(gPtr + 1966080) = rPtr[30];
            rPtr_4[30].x += rPtr[30].x;
            rPtr_4[30].y += rPtr[30].y;
            
            *(gPtr + 2031616) = rPtr[31];
            rPtr_4[31].x += rPtr[31].x;
            rPtr_4[31].y += rPtr[31].y;
            
        if(bid_cnt==thread_bs)
        
        {
        
        // 1's vector
        // tmp.x = (tx / 4 == 0) ? (rPtr_3[0].y + rPtr_3[0].x) * 2048: 0;
        // tmp.y = (tx / 4 == 0) ? (abs(rPtr_3[0].y) + abs(rPtr_3[0].x)) * 2048: 0;
        tmp = tmp_1;
        tmp_1.y += tmp.x;
        tmp_1.x = (abs(tmp.y) + abs(tmp.x));
        
        // 1's vector
        // tmp.x = (tx / 4 == 0) ? tmp_3.x : 0;
        tmp.x = tmp_3.x;
        tmp_3.y = tmp.x + tmp_3.y;
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 16, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 8, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 4, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 2, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 1, 32);
        
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 16, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 8, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 4, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 2, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 1, 32);

         // ToDo: can be optimized __shfl_sync
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 16, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 8, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 4, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 2, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 1, 32);
        __syncthreads();
        shPtr[(tx / 32) * 2] = tmp_1;
        shPtr[(tx / 32) * 2 + 1] = tmp_3;
        __syncthreads();
        
            tmp_1 = shPtr[(tx % 8) * 2];
            tmp_3 = shPtr[(tx % 8) * 2 + 1];
        
                tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 4, 32);
                tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 4, 32);
                tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 4, 32);
        
                tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 2, 32);
                tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 2, 32);
                tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 2, 32);
        
                tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 1, 32);
                tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 1, 32);
                tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 1, 32);
        
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(tx == 0)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if((blockIdx.x % thread_bs + 1) != round(abs(tmp_3.y) / abs(tmp_1.y)) && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 )  printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            //                                         bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x, tmp_3.y, tmp_3.y / tmp_1.y);
            // if(abs(tmp_1.y / tmp_1.x) > 0.001)printf("1, bid=%d bx=%d, by=%d, tx=%d: %f, %f, %f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // k = abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 ? bid : k;
            k = abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 ? round(abs(tmp_3.y) / abs(tmp_1.y)) : k;
            // k = abs(tmp_1.y) > 10 ? bid : k;
            // if(tx == 0) *(gPtr) = tmp_1;
            // if(tx == 0 && abs(tmp_1.y / tmp_1.x) > 1e-3)
            
            }
            // }            
            
    }
    
}

#include "../../../TurboFFT_radix_2_template.h"
template<>
__global__ void fft_radix_2<float2, 21, 1, 1, 1>(float2* inputs, float2* outputs, float2* twiddle, float2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    float2* shared = (float2*) ext_shared;
    int threadblock_per_SM = 4;
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
    
    rPtr_2[0] = *(checksum_DFT + 2048 - 2 + tx + 0);
    shPtr[tx + 0] = rPtr_2[0];
    
    rPtr_2[1] = *(checksum_DFT + 2048 - 2 + tx + 256);
    shPtr[tx + 256] = rPtr_2[1];
    
    rPtr_2[2] = *(checksum_DFT + 2048 - 2 + tx + 512);
    shPtr[tx + 512] = rPtr_2[2];
    
    rPtr_2[3] = *(checksum_DFT + 2048 - 2 + tx + 768);
    shPtr[tx + 768] = rPtr_2[3];
    
    rPtr_2[4] = *(checksum_DFT + 2048 - 2 + tx + 1024);
    shPtr[tx + 1024] = rPtr_2[4];
    
    rPtr_2[5] = *(checksum_DFT + 2048 - 2 + tx + 1280);
    shPtr[tx + 1280] = rPtr_2[5];
    
    rPtr_2[6] = *(checksum_DFT + 2048 - 2 + tx + 1536);
    shPtr[tx + 1536] = rPtr_2[6];
    
    rPtr_2[7] = *(checksum_DFT + 2048 - 2 + tx + 1792);
    shPtr[tx + 1792] = rPtr_2[7];
    
    __syncthreads();
    tmp_1.x = 0;
    tmp_1.y = 0;
    tmp_2.x = 0;
    tmp_2.y = 0;
    tmp_3.x = 0;
    tmp_3.y = 0;
    
    rPtr_2[0] = *(shPtr +  tx / 4 + 0);
    rPtr_3[0].x = 0; rPtr_3[0].y = 0;
    rPtr_4[0].x = 0; rPtr_4[0].y = 0;
    
    rPtr_2[1] = *(shPtr +  tx / 4 + 64);
    rPtr_3[1].x = 0; rPtr_3[1].y = 0;
    rPtr_4[1].x = 0; rPtr_4[1].y = 0;
    
    rPtr_2[2] = *(shPtr +  tx / 4 + 128);
    rPtr_3[2].x = 0; rPtr_3[2].y = 0;
    rPtr_4[2].x = 0; rPtr_4[2].y = 0;
    
    rPtr_2[3] = *(shPtr +  tx / 4 + 192);
    rPtr_3[3].x = 0; rPtr_3[3].y = 0;
    rPtr_4[3].x = 0; rPtr_4[3].y = 0;
    
    rPtr_2[4] = *(shPtr +  tx / 4 + 256);
    rPtr_3[4].x = 0; rPtr_3[4].y = 0;
    rPtr_4[4].x = 0; rPtr_4[4].y = 0;
    
    rPtr_2[5] = *(shPtr +  tx / 4 + 320);
    rPtr_3[5].x = 0; rPtr_3[5].y = 0;
    rPtr_4[5].x = 0; rPtr_4[5].y = 0;
    
    rPtr_2[6] = *(shPtr +  tx / 4 + 384);
    rPtr_3[6].x = 0; rPtr_3[6].y = 0;
    rPtr_4[6].x = 0; rPtr_4[6].y = 0;
    
    rPtr_2[7] = *(shPtr +  tx / 4 + 448);
    rPtr_3[7].x = 0; rPtr_3[7].y = 0;
    rPtr_4[7].x = 0; rPtr_4[7].y = 0;
    
    rPtr_2[8] = *(shPtr +  tx / 4 + 512);
    rPtr_3[8].x = 0; rPtr_3[8].y = 0;
    rPtr_4[8].x = 0; rPtr_4[8].y = 0;
    
    rPtr_2[9] = *(shPtr +  tx / 4 + 576);
    rPtr_3[9].x = 0; rPtr_3[9].y = 0;
    rPtr_4[9].x = 0; rPtr_4[9].y = 0;
    
    rPtr_2[10] = *(shPtr +  tx / 4 + 640);
    rPtr_3[10].x = 0; rPtr_3[10].y = 0;
    rPtr_4[10].x = 0; rPtr_4[10].y = 0;
    
    rPtr_2[11] = *(shPtr +  tx / 4 + 704);
    rPtr_3[11].x = 0; rPtr_3[11].y = 0;
    rPtr_4[11].x = 0; rPtr_4[11].y = 0;
    
    rPtr_2[12] = *(shPtr +  tx / 4 + 768);
    rPtr_3[12].x = 0; rPtr_3[12].y = 0;
    rPtr_4[12].x = 0; rPtr_4[12].y = 0;
    
    rPtr_2[13] = *(shPtr +  tx / 4 + 832);
    rPtr_3[13].x = 0; rPtr_3[13].y = 0;
    rPtr_4[13].x = 0; rPtr_4[13].y = 0;
    
    rPtr_2[14] = *(shPtr +  tx / 4 + 896);
    rPtr_3[14].x = 0; rPtr_3[14].y = 0;
    rPtr_4[14].x = 0; rPtr_4[14].y = 0;
    
    rPtr_2[15] = *(shPtr +  tx / 4 + 960);
    rPtr_3[15].x = 0; rPtr_3[15].y = 0;
    rPtr_4[15].x = 0; rPtr_4[15].y = 0;
    
    rPtr_2[16] = *(shPtr +  tx / 4 + 1024);
    rPtr_3[16].x = 0; rPtr_3[16].y = 0;
    rPtr_4[16].x = 0; rPtr_4[16].y = 0;
    
    rPtr_2[17] = *(shPtr +  tx / 4 + 1088);
    rPtr_3[17].x = 0; rPtr_3[17].y = 0;
    rPtr_4[17].x = 0; rPtr_4[17].y = 0;
    
    rPtr_2[18] = *(shPtr +  tx / 4 + 1152);
    rPtr_3[18].x = 0; rPtr_3[18].y = 0;
    rPtr_4[18].x = 0; rPtr_4[18].y = 0;
    
    rPtr_2[19] = *(shPtr +  tx / 4 + 1216);
    rPtr_3[19].x = 0; rPtr_3[19].y = 0;
    rPtr_4[19].x = 0; rPtr_4[19].y = 0;
    
    rPtr_2[20] = *(shPtr +  tx / 4 + 1280);
    rPtr_3[20].x = 0; rPtr_3[20].y = 0;
    rPtr_4[20].x = 0; rPtr_4[20].y = 0;
    
    rPtr_2[21] = *(shPtr +  tx / 4 + 1344);
    rPtr_3[21].x = 0; rPtr_3[21].y = 0;
    rPtr_4[21].x = 0; rPtr_4[21].y = 0;
    
    rPtr_2[22] = *(shPtr +  tx / 4 + 1408);
    rPtr_3[22].x = 0; rPtr_3[22].y = 0;
    rPtr_4[22].x = 0; rPtr_4[22].y = 0;
    
    rPtr_2[23] = *(shPtr +  tx / 4 + 1472);
    rPtr_3[23].x = 0; rPtr_3[23].y = 0;
    rPtr_4[23].x = 0; rPtr_4[23].y = 0;
    
    rPtr_2[24] = *(shPtr +  tx / 4 + 1536);
    rPtr_3[24].x = 0; rPtr_3[24].y = 0;
    rPtr_4[24].x = 0; rPtr_4[24].y = 0;
    
    rPtr_2[25] = *(shPtr +  tx / 4 + 1600);
    rPtr_3[25].x = 0; rPtr_3[25].y = 0;
    rPtr_4[25].x = 0; rPtr_4[25].y = 0;
    
    rPtr_2[26] = *(shPtr +  tx / 4 + 1664);
    rPtr_3[26].x = 0; rPtr_3[26].y = 0;
    rPtr_4[26].x = 0; rPtr_4[26].y = 0;
    
    rPtr_2[27] = *(shPtr +  tx / 4 + 1728);
    rPtr_3[27].x = 0; rPtr_3[27].y = 0;
    rPtr_4[27].x = 0; rPtr_4[27].y = 0;
    
    rPtr_2[28] = *(shPtr +  tx / 4 + 1792);
    rPtr_3[28].x = 0; rPtr_3[28].y = 0;
    rPtr_4[28].x = 0; rPtr_4[28].y = 0;
    
    rPtr_2[29] = *(shPtr +  tx / 4 + 1856);
    rPtr_3[29].x = 0; rPtr_3[29].y = 0;
    rPtr_4[29].x = 0; rPtr_4[29].y = 0;
    
    rPtr_2[30] = *(shPtr +  tx / 4 + 1920);
    rPtr_3[30].x = 0; rPtr_3[30].y = 0;
    rPtr_4[30].x = 0; rPtr_4[30].y = 0;
    
    rPtr_2[31] = *(shPtr +  tx / 4 + 1984);
    rPtr_3[31].x = 0; rPtr_3[31].y = 0;
    rPtr_4[31].x = 0; rPtr_4[31].y = 0;
    
    __syncthreads();
    int bid = 0;
    for(bid = (blockIdx.x / tb_gap) * tb_gap * thread_bs + blockIdx.x % tb_gap;
                bid_cnt < thread_bs && bid < (2097152 * BS + 8192 - 1) / 8192; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 4 * 1;
    
    gPtr += (bx % 1) * 2048 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 2048;
    bx = bx / 256;
    
    gPtr += tx % 4 * 2048;
    
    gPtr += (bx % BS * 2097152);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        // tmp = checksum_DFT[tx / 4 + 0];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[0], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[0], rPtr_2[0])
        turboFFT_ZMUL(tmp, rPtr[0], rPtr_2[0])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[1] = *(gPtr + 64);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
        // tmp = checksum_DFT[tx / 4 + 64];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[1], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[1], rPtr_2[1])
        turboFFT_ZMUL(tmp, rPtr[1], rPtr_2[1])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[2] = *(gPtr + 128);
        rPtr_3[2].x += rPtr[2].x;
        rPtr_3[2].y += rPtr[2].y;
        
        // tmp = checksum_DFT[tx / 4 + 128];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[2], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[2], rPtr_2[2])
        turboFFT_ZMUL(tmp, rPtr[2], rPtr_2[2])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[3] = *(gPtr + 192);
        rPtr_3[3].x += rPtr[3].x;
        rPtr_3[3].y += rPtr[3].y;
        
        // tmp = checksum_DFT[tx / 4 + 192];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[3], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[3], rPtr_2[3])
        turboFFT_ZMUL(tmp, rPtr[3], rPtr_2[3])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[4] = *(gPtr + 256);
        rPtr_3[4].x += rPtr[4].x;
        rPtr_3[4].y += rPtr[4].y;
        
        // tmp = checksum_DFT[tx / 4 + 256];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[4], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[4], rPtr_2[4])
        turboFFT_ZMUL(tmp, rPtr[4], rPtr_2[4])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[5] = *(gPtr + 320);
        rPtr_3[5].x += rPtr[5].x;
        rPtr_3[5].y += rPtr[5].y;
        
        // tmp = checksum_DFT[tx / 4 + 320];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[5], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[5], rPtr_2[5])
        turboFFT_ZMUL(tmp, rPtr[5], rPtr_2[5])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[6] = *(gPtr + 384);
        rPtr_3[6].x += rPtr[6].x;
        rPtr_3[6].y += rPtr[6].y;
        
        // tmp = checksum_DFT[tx / 4 + 384];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[6], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[6], rPtr_2[6])
        turboFFT_ZMUL(tmp, rPtr[6], rPtr_2[6])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[7] = *(gPtr + 448);
        rPtr_3[7].x += rPtr[7].x;
        rPtr_3[7].y += rPtr[7].y;
        
        // tmp = checksum_DFT[tx / 4 + 448];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[7], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[7], rPtr_2[7])
        turboFFT_ZMUL(tmp, rPtr[7], rPtr_2[7])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[8] = *(gPtr + 512);
        rPtr_3[8].x += rPtr[8].x;
        rPtr_3[8].y += rPtr[8].y;
        
        // tmp = checksum_DFT[tx / 4 + 512];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[8], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[8], rPtr_2[8])
        turboFFT_ZMUL(tmp, rPtr[8], rPtr_2[8])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[9] = *(gPtr + 576);
        rPtr_3[9].x += rPtr[9].x;
        rPtr_3[9].y += rPtr[9].y;
        
        // tmp = checksum_DFT[tx / 4 + 576];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[9], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[9], rPtr_2[9])
        turboFFT_ZMUL(tmp, rPtr[9], rPtr_2[9])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[10] = *(gPtr + 640);
        rPtr_3[10].x += rPtr[10].x;
        rPtr_3[10].y += rPtr[10].y;
        
        // tmp = checksum_DFT[tx / 4 + 640];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[10], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[10], rPtr_2[10])
        turboFFT_ZMUL(tmp, rPtr[10], rPtr_2[10])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[11] = *(gPtr + 704);
        rPtr_3[11].x += rPtr[11].x;
        rPtr_3[11].y += rPtr[11].y;
        
        // tmp = checksum_DFT[tx / 4 + 704];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[11], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[11], rPtr_2[11])
        turboFFT_ZMUL(tmp, rPtr[11], rPtr_2[11])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[12] = *(gPtr + 768);
        rPtr_3[12].x += rPtr[12].x;
        rPtr_3[12].y += rPtr[12].y;
        
        // tmp = checksum_DFT[tx / 4 + 768];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[12], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[12], rPtr_2[12])
        turboFFT_ZMUL(tmp, rPtr[12], rPtr_2[12])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[13] = *(gPtr + 832);
        rPtr_3[13].x += rPtr[13].x;
        rPtr_3[13].y += rPtr[13].y;
        
        // tmp = checksum_DFT[tx / 4 + 832];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[13], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[13], rPtr_2[13])
        turboFFT_ZMUL(tmp, rPtr[13], rPtr_2[13])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[14] = *(gPtr + 896);
        rPtr_3[14].x += rPtr[14].x;
        rPtr_3[14].y += rPtr[14].y;
        
        // tmp = checksum_DFT[tx / 4 + 896];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[14], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[14], rPtr_2[14])
        turboFFT_ZMUL(tmp, rPtr[14], rPtr_2[14])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[15] = *(gPtr + 960);
        rPtr_3[15].x += rPtr[15].x;
        rPtr_3[15].y += rPtr[15].y;
        
        // tmp = checksum_DFT[tx / 4 + 960];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[15], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[15], rPtr_2[15])
        turboFFT_ZMUL(tmp, rPtr[15], rPtr_2[15])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[16] = *(gPtr + 1024);
        rPtr_3[16].x += rPtr[16].x;
        rPtr_3[16].y += rPtr[16].y;
        
        // tmp = checksum_DFT[tx / 4 + 1024];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[16], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[16], rPtr_2[16])
        turboFFT_ZMUL(tmp, rPtr[16], rPtr_2[16])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[17] = *(gPtr + 1088);
        rPtr_3[17].x += rPtr[17].x;
        rPtr_3[17].y += rPtr[17].y;
        
        // tmp = checksum_DFT[tx / 4 + 1088];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[17], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[17], rPtr_2[17])
        turboFFT_ZMUL(tmp, rPtr[17], rPtr_2[17])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[18] = *(gPtr + 1152);
        rPtr_3[18].x += rPtr[18].x;
        rPtr_3[18].y += rPtr[18].y;
        
        // tmp = checksum_DFT[tx / 4 + 1152];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[18], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[18], rPtr_2[18])
        turboFFT_ZMUL(tmp, rPtr[18], rPtr_2[18])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[19] = *(gPtr + 1216);
        rPtr_3[19].x += rPtr[19].x;
        rPtr_3[19].y += rPtr[19].y;
        
        // tmp = checksum_DFT[tx / 4 + 1216];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[19], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[19], rPtr_2[19])
        turboFFT_ZMUL(tmp, rPtr[19], rPtr_2[19])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[20] = *(gPtr + 1280);
        rPtr_3[20].x += rPtr[20].x;
        rPtr_3[20].y += rPtr[20].y;
        
        // tmp = checksum_DFT[tx / 4 + 1280];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[20], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[20], rPtr_2[20])
        turboFFT_ZMUL(tmp, rPtr[20], rPtr_2[20])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[21] = *(gPtr + 1344);
        rPtr_3[21].x += rPtr[21].x;
        rPtr_3[21].y += rPtr[21].y;
        
        // tmp = checksum_DFT[tx / 4 + 1344];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[21], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[21], rPtr_2[21])
        turboFFT_ZMUL(tmp, rPtr[21], rPtr_2[21])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[22] = *(gPtr + 1408);
        rPtr_3[22].x += rPtr[22].x;
        rPtr_3[22].y += rPtr[22].y;
        
        // tmp = checksum_DFT[tx / 4 + 1408];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[22], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[22], rPtr_2[22])
        turboFFT_ZMUL(tmp, rPtr[22], rPtr_2[22])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[23] = *(gPtr + 1472);
        rPtr_3[23].x += rPtr[23].x;
        rPtr_3[23].y += rPtr[23].y;
        
        // tmp = checksum_DFT[tx / 4 + 1472];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[23], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[23], rPtr_2[23])
        turboFFT_ZMUL(tmp, rPtr[23], rPtr_2[23])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[24] = *(gPtr + 1536);
        rPtr_3[24].x += rPtr[24].x;
        rPtr_3[24].y += rPtr[24].y;
        
        // tmp = checksum_DFT[tx / 4 + 1536];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[24], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[24], rPtr_2[24])
        turboFFT_ZMUL(tmp, rPtr[24], rPtr_2[24])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[25] = *(gPtr + 1600);
        rPtr_3[25].x += rPtr[25].x;
        rPtr_3[25].y += rPtr[25].y;
        
        // tmp = checksum_DFT[tx / 4 + 1600];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[25], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[25], rPtr_2[25])
        turboFFT_ZMUL(tmp, rPtr[25], rPtr_2[25])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[26] = *(gPtr + 1664);
        rPtr_3[26].x += rPtr[26].x;
        rPtr_3[26].y += rPtr[26].y;
        
        // tmp = checksum_DFT[tx / 4 + 1664];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[26], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[26], rPtr_2[26])
        turboFFT_ZMUL(tmp, rPtr[26], rPtr_2[26])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[27] = *(gPtr + 1728);
        rPtr_3[27].x += rPtr[27].x;
        rPtr_3[27].y += rPtr[27].y;
        
        // tmp = checksum_DFT[tx / 4 + 1728];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[27], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[27], rPtr_2[27])
        turboFFT_ZMUL(tmp, rPtr[27], rPtr_2[27])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[28] = *(gPtr + 1792);
        rPtr_3[28].x += rPtr[28].x;
        rPtr_3[28].y += rPtr[28].y;
        
        // tmp = checksum_DFT[tx / 4 + 1792];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[28], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[28], rPtr_2[28])
        turboFFT_ZMUL(tmp, rPtr[28], rPtr_2[28])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[29] = *(gPtr + 1856);
        rPtr_3[29].x += rPtr[29].x;
        rPtr_3[29].y += rPtr[29].y;
        
        // tmp = checksum_DFT[tx / 4 + 1856];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[29], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[29], rPtr_2[29])
        turboFFT_ZMUL(tmp, rPtr[29], rPtr_2[29])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[30] = *(gPtr + 1920);
        rPtr_3[30].x += rPtr[30].x;
        rPtr_3[30].y += rPtr[30].y;
        
        // tmp = checksum_DFT[tx / 4 + 1920];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[30], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[30], rPtr_2[30])
        turboFFT_ZMUL(tmp, rPtr[30], rPtr_2[30])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[31] = *(gPtr + 1984);
        rPtr_3[31].x += rPtr[31].x;
        rPtr_3[31].y += rPtr[31].y;
        
        // tmp = checksum_DFT[tx / 4 + 1984];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[31], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[31], rPtr_2[31])
        turboFFT_ZMUL(tmp, rPtr[31], rPtr_2[31])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        // tmp_3.x += bid_cnt * (rPtr[0].x + rPtr[0].y) * 2048;
        
        rPtr[0].x += (threadIdx.x == 0 && bid_cnt == (blockIdx.x % thread_bs + 1)) ? 100: 0;
        
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
    
    offset += ((tx / 1) % 4) * 1;
    
    j = tx / 4;
    
    offset += ((tx / 4) % 2) * 128;
    
    offset += ((tx / 8) % 32) * 256;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.003067961661145091f);
    delta_angle.y = __sinf(j * -0.003067961661145091f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 4] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 8] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 12] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 16] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 20] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 24] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 28] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 32] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 36] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 40] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 44] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 48] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 52] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 56] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 60] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 64] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 68] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 72] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 76] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 80] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 84] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 88] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 92] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 96] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 100] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 104] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 108] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 112] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 116] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 120] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 124] = rPtr[31];
    
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
    
    j = 0;
    offset  = 0;
    
    offset += ((tx / 1) % 4) * 1;
    
    offset += ((tx / 4) % 32) * 4;
    
    j = tx / 128;
    
    offset += ((tx / 128) % 2) * 4096;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.09817477315664291f);
    delta_angle.y = __sinf(j * -0.09817477315664291f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 128] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 256] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 384] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 512] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 640] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 768] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 896] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 1024] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 1152] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 1280] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 1408] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 1536] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 1664] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 1792] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 1920] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 2048] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 2176] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 2304] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 2432] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 2560] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 2688] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 2816] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 2944] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 3072] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 3200] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 3328] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 3456] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 3584] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 3712] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 3840] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 3968] = rPtr[31];
    
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
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
            
    bx = bid;
    tx = threadIdx.x;
    gPtr = outputs;
    
    gPtr += tx / 4 * 1024;
    
    gPtr += (bx % 1) * 2048 * 1024;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 1;
    bx = bx / 256;
    
    gPtr += tx % 4 * 1;
    
    gPtr += (bx % BS * 2097152);
    
        // 1's vector
        // tmp_3.y -=  (rPtr[0].y + rPtr[0].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[0].y + rPtr[0].x);
        turboFFT_ZMUL(tmp, rPtr[0],r[(0 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[0], r[(0 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[0], r[(0 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[1].y + rPtr[1].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[1].y + rPtr[1].x);
        turboFFT_ZMUL(tmp, rPtr[1],r[(64 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[1], r[(64 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[1], r[(64 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[2].y + rPtr[2].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[2].y + rPtr[2].x);
        turboFFT_ZMUL(tmp, rPtr[2],r[(128 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[2], r[(128 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[2], r[(128 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[3].y + rPtr[3].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[3].y + rPtr[3].x);
        turboFFT_ZMUL(tmp, rPtr[3],r[(192 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[3], r[(192 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[3], r[(192 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[4].y + rPtr[4].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[4].y + rPtr[4].x);
        turboFFT_ZMUL(tmp, rPtr[4],r[(256 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[4], r[(256 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[4], r[(256 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[5].y + rPtr[5].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[5].y + rPtr[5].x);
        turboFFT_ZMUL(tmp, rPtr[5],r[(320 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[5], r[(320 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[5], r[(320 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[6].y + rPtr[6].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[6].y + rPtr[6].x);
        turboFFT_ZMUL(tmp, rPtr[6],r[(384 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[6], r[(384 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[6], r[(384 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[7].y + rPtr[7].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[7].y + rPtr[7].x);
        turboFFT_ZMUL(tmp, rPtr[7],r[(448 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[7], r[(448 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[7], r[(448 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[8].y + rPtr[8].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[8].y + rPtr[8].x);
        turboFFT_ZMUL(tmp, rPtr[8],r[(512 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[8], r[(512 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[8], r[(512 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[9].y + rPtr[9].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[9].y + rPtr[9].x);
        turboFFT_ZMUL(tmp, rPtr[9],r[(576 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[9], r[(576 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[9], r[(576 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[10].y + rPtr[10].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[10].y + rPtr[10].x);
        turboFFT_ZMUL(tmp, rPtr[10],r[(640 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[10], r[(640 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[10], r[(640 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[11].y + rPtr[11].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[11].y + rPtr[11].x);
        turboFFT_ZMUL(tmp, rPtr[11],r[(704 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[11], r[(704 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[11], r[(704 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[12].y + rPtr[12].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[12].y + rPtr[12].x);
        turboFFT_ZMUL(tmp, rPtr[12],r[(768 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[12], r[(768 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[12], r[(768 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[13].y + rPtr[13].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[13].y + rPtr[13].x);
        turboFFT_ZMUL(tmp, rPtr[13],r[(832 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[13], r[(832 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[13], r[(832 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[14].y + rPtr[14].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[14].y + rPtr[14].x);
        turboFFT_ZMUL(tmp, rPtr[14],r[(896 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[14], r[(896 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[14], r[(896 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[15].y + rPtr[15].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[15].y + rPtr[15].x);
        turboFFT_ZMUL(tmp, rPtr[15],r[(960 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[15], r[(960 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[15], r[(960 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[16].y + rPtr[16].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[16].y + rPtr[16].x);
        turboFFT_ZMUL(tmp, rPtr[16],r[(1024 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[16], r[(1024 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[16], r[(1024 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[17].y + rPtr[17].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[17].y + rPtr[17].x);
        turboFFT_ZMUL(tmp, rPtr[17],r[(1088 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[17], r[(1088 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[17], r[(1088 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[18].y + rPtr[18].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[18].y + rPtr[18].x);
        turboFFT_ZMUL(tmp, rPtr[18],r[(1152 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[18], r[(1152 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[18], r[(1152 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[19].y + rPtr[19].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[19].y + rPtr[19].x);
        turboFFT_ZMUL(tmp, rPtr[19],r[(1216 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[19], r[(1216 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[19], r[(1216 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[20].y + rPtr[20].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[20].y + rPtr[20].x);
        turboFFT_ZMUL(tmp, rPtr[20],r[(1280 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[20], r[(1280 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[20], r[(1280 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[21].y + rPtr[21].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[21].y + rPtr[21].x);
        turboFFT_ZMUL(tmp, rPtr[21],r[(1344 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[21], r[(1344 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[21], r[(1344 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[22].y + rPtr[22].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[22].y + rPtr[22].x);
        turboFFT_ZMUL(tmp, rPtr[22],r[(1408 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[22], r[(1408 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[22], r[(1408 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[23].y + rPtr[23].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[23].y + rPtr[23].x);
        turboFFT_ZMUL(tmp, rPtr[23],r[(1472 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[23], r[(1472 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[23], r[(1472 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[24].y + rPtr[24].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[24].y + rPtr[24].x);
        turboFFT_ZMUL(tmp, rPtr[24],r[(1536 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[24], r[(1536 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[24], r[(1536 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[25].y + rPtr[25].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[25].y + rPtr[25].x);
        turboFFT_ZMUL(tmp, rPtr[25],r[(1600 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[25], r[(1600 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[25], r[(1600 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[26].y + rPtr[26].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[26].y + rPtr[26].x);
        turboFFT_ZMUL(tmp, rPtr[26],r[(1664 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[26], r[(1664 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[26], r[(1664 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[27].y + rPtr[27].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[27].y + rPtr[27].x);
        turboFFT_ZMUL(tmp, rPtr[27],r[(1728 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[27], r[(1728 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[27], r[(1728 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[28].y + rPtr[28].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[28].y + rPtr[28].x);
        turboFFT_ZMUL(tmp, rPtr[28],r[(1792 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[28], r[(1792 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[28], r[(1792 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[29].y + rPtr[29].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[29].y + rPtr[29].x);
        turboFFT_ZMUL(tmp, rPtr[29],r[(1856 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[29], r[(1856 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[29], r[(1856 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[30].y + rPtr[30].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[30].y + rPtr[30].x);
        turboFFT_ZMUL(tmp, rPtr[30],r[(1920 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[30], r[(1920 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[30], r[(1920 + tx / 4) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[31].y + rPtr[31].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[31].y + rPtr[31].x);
        turboFFT_ZMUL(tmp, rPtr[31],r[(1984 + tx / 4) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[31], r[(1984 + tx / 4) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[31], r[(1984 + tx / 4) % 3])
        
            *(gPtr + 0) = rPtr[0];
            rPtr_4[0].x += rPtr[0].x;
            rPtr_4[0].y += rPtr[0].y;
            
            *(gPtr + 65536) = rPtr[1];
            rPtr_4[1].x += rPtr[1].x;
            rPtr_4[1].y += rPtr[1].y;
            
            *(gPtr + 131072) = rPtr[2];
            rPtr_4[2].x += rPtr[2].x;
            rPtr_4[2].y += rPtr[2].y;
            
            *(gPtr + 196608) = rPtr[3];
            rPtr_4[3].x += rPtr[3].x;
            rPtr_4[3].y += rPtr[3].y;
            
            *(gPtr + 262144) = rPtr[4];
            rPtr_4[4].x += rPtr[4].x;
            rPtr_4[4].y += rPtr[4].y;
            
            *(gPtr + 327680) = rPtr[5];
            rPtr_4[5].x += rPtr[5].x;
            rPtr_4[5].y += rPtr[5].y;
            
            *(gPtr + 393216) = rPtr[6];
            rPtr_4[6].x += rPtr[6].x;
            rPtr_4[6].y += rPtr[6].y;
            
            *(gPtr + 458752) = rPtr[7];
            rPtr_4[7].x += rPtr[7].x;
            rPtr_4[7].y += rPtr[7].y;
            
            *(gPtr + 524288) = rPtr[8];
            rPtr_4[8].x += rPtr[8].x;
            rPtr_4[8].y += rPtr[8].y;
            
            *(gPtr + 589824) = rPtr[9];
            rPtr_4[9].x += rPtr[9].x;
            rPtr_4[9].y += rPtr[9].y;
            
            *(gPtr + 655360) = rPtr[10];
            rPtr_4[10].x += rPtr[10].x;
            rPtr_4[10].y += rPtr[10].y;
            
            *(gPtr + 720896) = rPtr[11];
            rPtr_4[11].x += rPtr[11].x;
            rPtr_4[11].y += rPtr[11].y;
            
            *(gPtr + 786432) = rPtr[12];
            rPtr_4[12].x += rPtr[12].x;
            rPtr_4[12].y += rPtr[12].y;
            
            *(gPtr + 851968) = rPtr[13];
            rPtr_4[13].x += rPtr[13].x;
            rPtr_4[13].y += rPtr[13].y;
            
            *(gPtr + 917504) = rPtr[14];
            rPtr_4[14].x += rPtr[14].x;
            rPtr_4[14].y += rPtr[14].y;
            
            *(gPtr + 983040) = rPtr[15];
            rPtr_4[15].x += rPtr[15].x;
            rPtr_4[15].y += rPtr[15].y;
            
            *(gPtr + 1048576) = rPtr[16];
            rPtr_4[16].x += rPtr[16].x;
            rPtr_4[16].y += rPtr[16].y;
            
            *(gPtr + 1114112) = rPtr[17];
            rPtr_4[17].x += rPtr[17].x;
            rPtr_4[17].y += rPtr[17].y;
            
            *(gPtr + 1179648) = rPtr[18];
            rPtr_4[18].x += rPtr[18].x;
            rPtr_4[18].y += rPtr[18].y;
            
            *(gPtr + 1245184) = rPtr[19];
            rPtr_4[19].x += rPtr[19].x;
            rPtr_4[19].y += rPtr[19].y;
            
            *(gPtr + 1310720) = rPtr[20];
            rPtr_4[20].x += rPtr[20].x;
            rPtr_4[20].y += rPtr[20].y;
            
            *(gPtr + 1376256) = rPtr[21];
            rPtr_4[21].x += rPtr[21].x;
            rPtr_4[21].y += rPtr[21].y;
            
            *(gPtr + 1441792) = rPtr[22];
            rPtr_4[22].x += rPtr[22].x;
            rPtr_4[22].y += rPtr[22].y;
            
            *(gPtr + 1507328) = rPtr[23];
            rPtr_4[23].x += rPtr[23].x;
            rPtr_4[23].y += rPtr[23].y;
            
            *(gPtr + 1572864) = rPtr[24];
            rPtr_4[24].x += rPtr[24].x;
            rPtr_4[24].y += rPtr[24].y;
            
            *(gPtr + 1638400) = rPtr[25];
            rPtr_4[25].x += rPtr[25].x;
            rPtr_4[25].y += rPtr[25].y;
            
            *(gPtr + 1703936) = rPtr[26];
            rPtr_4[26].x += rPtr[26].x;
            rPtr_4[26].y += rPtr[26].y;
            
            *(gPtr + 1769472) = rPtr[27];
            rPtr_4[27].x += rPtr[27].x;
            rPtr_4[27].y += rPtr[27].y;
            
            *(gPtr + 1835008) = rPtr[28];
            rPtr_4[28].x += rPtr[28].x;
            rPtr_4[28].y += rPtr[28].y;
            
            *(gPtr + 1900544) = rPtr[29];
            rPtr_4[29].x += rPtr[29].x;
            rPtr_4[29].y += rPtr[29].y;
            
            *(gPtr + 1966080) = rPtr[30];
            rPtr_4[30].x += rPtr[30].x;
            rPtr_4[30].y += rPtr[30].y;
            
            *(gPtr + 2031616) = rPtr[31];
            rPtr_4[31].x += rPtr[31].x;
            rPtr_4[31].y += rPtr[31].y;
            
        if(bid_cnt==thread_bs)
        
        {
        
        // 1's vector
        // tmp.x = (tx / 4 == 0) ? (rPtr_3[0].y + rPtr_3[0].x) * 2048: 0;
        // tmp.y = (tx / 4 == 0) ? (abs(rPtr_3[0].y) + abs(rPtr_3[0].x)) * 2048: 0;
        tmp = tmp_1;
        tmp_1.y += tmp.x;
        tmp_1.x = (abs(tmp.y) + abs(tmp.x));
        
        // 1's vector
        // tmp.x = (tx / 4 == 0) ? tmp_3.x : 0;
        tmp.x = tmp_3.x;
        tmp_3.y = tmp.x + tmp_3.y;
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 16, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 8, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 4, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 2, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 1, 32);
        
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 16, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 8, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 4, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 2, 32);
        tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 1, 32);

         // ToDo: can be optimized __shfl_sync
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 16, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 8, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 4, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 2, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 1, 32);
        __syncthreads();
        shPtr[(tx / 32) * 2] = tmp_1;
        shPtr[(tx / 32) * 2 + 1] = tmp_3;
        __syncthreads();
        
            tmp_1 = shPtr[(tx % 8) * 2];
            tmp_3 = shPtr[(tx % 8) * 2 + 1];
        
                tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 4, 32);
                tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 4, 32);
                tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 4, 32);
        
                tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 2, 32);
                tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 2, 32);
                tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 2, 32);
        
                tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 1, 32);
                tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 1, 32);
                tmp_3.y += __shfl_xor_sync(0xffffffff, tmp_3.y, 1, 32);
        
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(tx == 0)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if((blockIdx.x % thread_bs + 1) != round(abs(tmp_3.y) / abs(tmp_1.y)) && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 )  printf("1, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            //                                         bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x, tmp_3.y, tmp_3.y / tmp_1.y);
            // if(abs(tmp_1.y / tmp_1.x) > 0.001)printf("1, bid=%d bx=%d, by=%d, tx=%d: %f, %f, %f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // k = abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 ? bid : k;
            k = abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 ? round(abs(tmp_3.y) / abs(tmp_1.y)) : k;
            // k = abs(tmp_1.y) > 10 ? bid : k;
            // if(tx == 0) *(gPtr) = tmp_1;
            // if(tx == 0 && abs(tmp_1.y / tmp_1.x) > 1e-3)
            
            }
            // }            
            
                }
                if(k != -1){
                
                bid = (blockIdx.x / tb_gap) * tb_gap * thread_bs + blockIdx.x % tb_gap + delta_bid * (k - 1);
                // if(threadIdx.x == 0)printf("bid=%d, upload=%d, bx=%d, tx=%d, k = %d\n", bid, 1, blockIdx.x, threadIdx.x, k);
                // bid = k;
                        
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 4 * 1;
    
    gPtr += (bx % 1) * 2048 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 2048;
    bx = bx / 256;
    
    gPtr += tx % 4 * 2048;
    
    gPtr += (bx % BS * 2097152);
    
        // rPtr[0] = rPtr_3[0];
        rPtr[0] = *(gPtr + 0);
        
        // rPtr[1] = rPtr_3[1];
        rPtr[1] = *(gPtr + 64);
        
        // rPtr[2] = rPtr_3[2];
        rPtr[2] = *(gPtr + 128);
        
        // rPtr[3] = rPtr_3[3];
        rPtr[3] = *(gPtr + 192);
        
        // rPtr[4] = rPtr_3[4];
        rPtr[4] = *(gPtr + 256);
        
        // rPtr[5] = rPtr_3[5];
        rPtr[5] = *(gPtr + 320);
        
        // rPtr[6] = rPtr_3[6];
        rPtr[6] = *(gPtr + 384);
        
        // rPtr[7] = rPtr_3[7];
        rPtr[7] = *(gPtr + 448);
        
        // rPtr[8] = rPtr_3[8];
        rPtr[8] = *(gPtr + 512);
        
        // rPtr[9] = rPtr_3[9];
        rPtr[9] = *(gPtr + 576);
        
        // rPtr[10] = rPtr_3[10];
        rPtr[10] = *(gPtr + 640);
        
        // rPtr[11] = rPtr_3[11];
        rPtr[11] = *(gPtr + 704);
        
        // rPtr[12] = rPtr_3[12];
        rPtr[12] = *(gPtr + 768);
        
        // rPtr[13] = rPtr_3[13];
        rPtr[13] = *(gPtr + 832);
        
        // rPtr[14] = rPtr_3[14];
        rPtr[14] = *(gPtr + 896);
        
        // rPtr[15] = rPtr_3[15];
        rPtr[15] = *(gPtr + 960);
        
        // rPtr[16] = rPtr_3[16];
        rPtr[16] = *(gPtr + 1024);
        
        // rPtr[17] = rPtr_3[17];
        rPtr[17] = *(gPtr + 1088);
        
        // rPtr[18] = rPtr_3[18];
        rPtr[18] = *(gPtr + 1152);
        
        // rPtr[19] = rPtr_3[19];
        rPtr[19] = *(gPtr + 1216);
        
        // rPtr[20] = rPtr_3[20];
        rPtr[20] = *(gPtr + 1280);
        
        // rPtr[21] = rPtr_3[21];
        rPtr[21] = *(gPtr + 1344);
        
        // rPtr[22] = rPtr_3[22];
        rPtr[22] = *(gPtr + 1408);
        
        // rPtr[23] = rPtr_3[23];
        rPtr[23] = *(gPtr + 1472);
        
        // rPtr[24] = rPtr_3[24];
        rPtr[24] = *(gPtr + 1536);
        
        // rPtr[25] = rPtr_3[25];
        rPtr[25] = *(gPtr + 1600);
        
        // rPtr[26] = rPtr_3[26];
        rPtr[26] = *(gPtr + 1664);
        
        // rPtr[27] = rPtr_3[27];
        rPtr[27] = *(gPtr + 1728);
        
        // rPtr[28] = rPtr_3[28];
        rPtr[28] = *(gPtr + 1792);
        
        // rPtr[29] = rPtr_3[29];
        rPtr[29] = *(gPtr + 1856);
        
        // rPtr[30] = rPtr_3[30];
        rPtr[30] = *(gPtr + 1920);
        
        // rPtr[31] = rPtr_3[31];
        rPtr[31] = *(gPtr + 1984);
        
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
    
    offset += ((tx / 1) % 4) * 1;
    
    j = tx / 4;
    
    offset += ((tx / 4) % 2) * 128;
    
    offset += ((tx / 8) % 32) * 256;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.003067961661145091f);
    delta_angle.y = __sinf(j * -0.003067961661145091f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 4] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 8] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 12] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 16] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 20] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 24] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 28] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 32] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 36] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 40] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 44] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 48] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 52] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 56] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 60] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 64] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 68] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 72] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 76] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 80] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 84] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 88] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 92] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 96] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 100] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 104] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 108] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 112] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 116] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 120] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 124] = rPtr[31];
    
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
    
    j = 0;
    offset  = 0;
    
    offset += ((tx / 1) % 4) * 1;
    
    offset += ((tx / 4) % 32) * 4;
    
    j = tx / 128;
    
    offset += ((tx / 128) % 2) * 4096;
    
    __syncthreads();
    
    delta_angle.x = __cosf(j * -0.09817477315664291f);
    delta_angle.y = __sinf(j * -0.09817477315664291f);
     
    angle.x = 1;
    angle.y = 0;
    
    shPtr[offset + 0] = rPtr[0];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[16];
    turboFFT_ZMUL(rPtr[16], tmp, angle);
    
    shPtr[offset + 128] = rPtr[16];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[8];
    turboFFT_ZMUL(rPtr[8], tmp, angle);
    
    shPtr[offset + 256] = rPtr[8];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[24];
    turboFFT_ZMUL(rPtr[24], tmp, angle);
    
    shPtr[offset + 384] = rPtr[24];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[4];
    turboFFT_ZMUL(rPtr[4], tmp, angle);
    
    shPtr[offset + 512] = rPtr[4];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[20];
    turboFFT_ZMUL(rPtr[20], tmp, angle);
    
    shPtr[offset + 640] = rPtr[20];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[12];
    turboFFT_ZMUL(rPtr[12], tmp, angle);
    
    shPtr[offset + 768] = rPtr[12];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[28];
    turboFFT_ZMUL(rPtr[28], tmp, angle);
    
    shPtr[offset + 896] = rPtr[28];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[2];
    turboFFT_ZMUL(rPtr[2], tmp, angle);
    
    shPtr[offset + 1024] = rPtr[2];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[18];
    turboFFT_ZMUL(rPtr[18], tmp, angle);
    
    shPtr[offset + 1152] = rPtr[18];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[10];
    turboFFT_ZMUL(rPtr[10], tmp, angle);
    
    shPtr[offset + 1280] = rPtr[10];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[26];
    turboFFT_ZMUL(rPtr[26], tmp, angle);
    
    shPtr[offset + 1408] = rPtr[26];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[6];
    turboFFT_ZMUL(rPtr[6], tmp, angle);
    
    shPtr[offset + 1536] = rPtr[6];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[22];
    turboFFT_ZMUL(rPtr[22], tmp, angle);
    
    shPtr[offset + 1664] = rPtr[22];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[14];
    turboFFT_ZMUL(rPtr[14], tmp, angle);
    
    shPtr[offset + 1792] = rPtr[14];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[30];
    turboFFT_ZMUL(rPtr[30], tmp, angle);
    
    shPtr[offset + 1920] = rPtr[30];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[1];
    turboFFT_ZMUL(rPtr[1], tmp, angle);
    
    shPtr[offset + 2048] = rPtr[1];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[17];
    turboFFT_ZMUL(rPtr[17], tmp, angle);
    
    shPtr[offset + 2176] = rPtr[17];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[9];
    turboFFT_ZMUL(rPtr[9], tmp, angle);
    
    shPtr[offset + 2304] = rPtr[9];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[25];
    turboFFT_ZMUL(rPtr[25], tmp, angle);
    
    shPtr[offset + 2432] = rPtr[25];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[5];
    turboFFT_ZMUL(rPtr[5], tmp, angle);
    
    shPtr[offset + 2560] = rPtr[5];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[21];
    turboFFT_ZMUL(rPtr[21], tmp, angle);
    
    shPtr[offset + 2688] = rPtr[21];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[13];
    turboFFT_ZMUL(rPtr[13], tmp, angle);
    
    shPtr[offset + 2816] = rPtr[13];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[29];
    turboFFT_ZMUL(rPtr[29], tmp, angle);
    
    shPtr[offset + 2944] = rPtr[29];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[3];
    turboFFT_ZMUL(rPtr[3], tmp, angle);
    
    shPtr[offset + 3072] = rPtr[3];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[19];
    turboFFT_ZMUL(rPtr[19], tmp, angle);
    
    shPtr[offset + 3200] = rPtr[19];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[11];
    turboFFT_ZMUL(rPtr[11], tmp, angle);
    
    shPtr[offset + 3328] = rPtr[11];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[27];
    turboFFT_ZMUL(rPtr[27], tmp, angle);
    
    shPtr[offset + 3456] = rPtr[27];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[7];
    turboFFT_ZMUL(rPtr[7], tmp, angle);
    
    shPtr[offset + 3584] = rPtr[7];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[23];
    turboFFT_ZMUL(rPtr[23], tmp, angle);
    
    shPtr[offset + 3712] = rPtr[23];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[15];
    turboFFT_ZMUL(rPtr[15], tmp, angle);
    
    shPtr[offset + 3840] = rPtr[15];
    
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = rPtr[31];
    turboFFT_ZMUL(rPtr[31], tmp, angle);
    
    shPtr[offset + 3968] = rPtr[31];
    
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
    
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
            
    bx = bid;
    tx = threadIdx.x;
    gPtr = outputs;
    
    gPtr += tx / 4 * 1024;
    
    gPtr += (bx % 1) * 2048 * 1024;
    bx = bx / 1;
    
    gPtr += (bx % 256) * 4 * 1;
    bx = bx / 256;
    
    gPtr += tx % 4 * 1;
    
    gPtr += (bx % BS * 2097152);
    
            // turboFFT_ZSUB(rPtr[0], rPtr[0], rPtr_4[0]);
            
            // turboFFT_ZSUB(rPtr[1], rPtr[1], rPtr_4[1]);
            
            // turboFFT_ZSUB(rPtr[2], rPtr[2], rPtr_4[2]);
            
            // turboFFT_ZSUB(rPtr[3], rPtr[3], rPtr_4[3]);
            
            // turboFFT_ZSUB(rPtr[4], rPtr[4], rPtr_4[4]);
            
            // turboFFT_ZSUB(rPtr[5], rPtr[5], rPtr_4[5]);
            
            // turboFFT_ZSUB(rPtr[6], rPtr[6], rPtr_4[6]);
            
            // turboFFT_ZSUB(rPtr[7], rPtr[7], rPtr_4[7]);
            
            // turboFFT_ZSUB(rPtr[8], rPtr[8], rPtr_4[8]);
            
            // turboFFT_ZSUB(rPtr[9], rPtr[9], rPtr_4[9]);
            
            // turboFFT_ZSUB(rPtr[10], rPtr[10], rPtr_4[10]);
            
            // turboFFT_ZSUB(rPtr[11], rPtr[11], rPtr_4[11]);
            
            // turboFFT_ZSUB(rPtr[12], rPtr[12], rPtr_4[12]);
            
            // turboFFT_ZSUB(rPtr[13], rPtr[13], rPtr_4[13]);
            
            // turboFFT_ZSUB(rPtr[14], rPtr[14], rPtr_4[14]);
            
            // turboFFT_ZSUB(rPtr[15], rPtr[15], rPtr_4[15]);
            
            // turboFFT_ZSUB(rPtr[16], rPtr[16], rPtr_4[16]);
            
            // turboFFT_ZSUB(rPtr[17], rPtr[17], rPtr_4[17]);
            
            // turboFFT_ZSUB(rPtr[18], rPtr[18], rPtr_4[18]);
            
            // turboFFT_ZSUB(rPtr[19], rPtr[19], rPtr_4[19]);
            
            // turboFFT_ZSUB(rPtr[20], rPtr[20], rPtr_4[20]);
            
            // turboFFT_ZSUB(rPtr[21], rPtr[21], rPtr_4[21]);
            
            // turboFFT_ZSUB(rPtr[22], rPtr[22], rPtr_4[22]);
            
            // turboFFT_ZSUB(rPtr[23], rPtr[23], rPtr_4[23]);
            
            // turboFFT_ZSUB(rPtr[24], rPtr[24], rPtr_4[24]);
            
            // turboFFT_ZSUB(rPtr[25], rPtr[25], rPtr_4[25]);
            
            // turboFFT_ZSUB(rPtr[26], rPtr[26], rPtr_4[26]);
            
            // turboFFT_ZSUB(rPtr[27], rPtr[27], rPtr_4[27]);
            
            // turboFFT_ZSUB(rPtr[28], rPtr[28], rPtr_4[28]);
            
            // turboFFT_ZSUB(rPtr[29], rPtr[29], rPtr_4[29]);
            
            // turboFFT_ZSUB(rPtr[30], rPtr[30], rPtr_4[30]);
            
            // turboFFT_ZSUB(rPtr[31], rPtr[31], rPtr_4[31]);
            
            // rPtr_3[0] = *(gPtr + 0);
            // turboFFT_ZADD(rPtr_3[0], rPtr_3[0], rPtr[0] );
            // *(gPtr + 0) = rPtr_3[0];
            *(gPtr + 0) = rPtr[0];
        
            // rPtr_3[1] = *(gPtr + 65536);
            // turboFFT_ZADD(rPtr_3[1], rPtr_3[1], rPtr[1] );
            // *(gPtr + 65536) = rPtr_3[1];
            *(gPtr + 65536) = rPtr[1];
        
            // rPtr_3[2] = *(gPtr + 131072);
            // turboFFT_ZADD(rPtr_3[2], rPtr_3[2], rPtr[2] );
            // *(gPtr + 131072) = rPtr_3[2];
            *(gPtr + 131072) = rPtr[2];
        
            // rPtr_3[3] = *(gPtr + 196608);
            // turboFFT_ZADD(rPtr_3[3], rPtr_3[3], rPtr[3] );
            // *(gPtr + 196608) = rPtr_3[3];
            *(gPtr + 196608) = rPtr[3];
        
            // rPtr_3[4] = *(gPtr + 262144);
            // turboFFT_ZADD(rPtr_3[4], rPtr_3[4], rPtr[4] );
            // *(gPtr + 262144) = rPtr_3[4];
            *(gPtr + 262144) = rPtr[4];
        
            // rPtr_3[5] = *(gPtr + 327680);
            // turboFFT_ZADD(rPtr_3[5], rPtr_3[5], rPtr[5] );
            // *(gPtr + 327680) = rPtr_3[5];
            *(gPtr + 327680) = rPtr[5];
        
            // rPtr_3[6] = *(gPtr + 393216);
            // turboFFT_ZADD(rPtr_3[6], rPtr_3[6], rPtr[6] );
            // *(gPtr + 393216) = rPtr_3[6];
            *(gPtr + 393216) = rPtr[6];
        
            // rPtr_3[7] = *(gPtr + 458752);
            // turboFFT_ZADD(rPtr_3[7], rPtr_3[7], rPtr[7] );
            // *(gPtr + 458752) = rPtr_3[7];
            *(gPtr + 458752) = rPtr[7];
        
            // rPtr_3[8] = *(gPtr + 524288);
            // turboFFT_ZADD(rPtr_3[8], rPtr_3[8], rPtr[8] );
            // *(gPtr + 524288) = rPtr_3[8];
            *(gPtr + 524288) = rPtr[8];
        
            // rPtr_3[9] = *(gPtr + 589824);
            // turboFFT_ZADD(rPtr_3[9], rPtr_3[9], rPtr[9] );
            // *(gPtr + 589824) = rPtr_3[9];
            *(gPtr + 589824) = rPtr[9];
        
            // rPtr_3[10] = *(gPtr + 655360);
            // turboFFT_ZADD(rPtr_3[10], rPtr_3[10], rPtr[10] );
            // *(gPtr + 655360) = rPtr_3[10];
            *(gPtr + 655360) = rPtr[10];
        
            // rPtr_3[11] = *(gPtr + 720896);
            // turboFFT_ZADD(rPtr_3[11], rPtr_3[11], rPtr[11] );
            // *(gPtr + 720896) = rPtr_3[11];
            *(gPtr + 720896) = rPtr[11];
        
            // rPtr_3[12] = *(gPtr + 786432);
            // turboFFT_ZADD(rPtr_3[12], rPtr_3[12], rPtr[12] );
            // *(gPtr + 786432) = rPtr_3[12];
            *(gPtr + 786432) = rPtr[12];
        
            // rPtr_3[13] = *(gPtr + 851968);
            // turboFFT_ZADD(rPtr_3[13], rPtr_3[13], rPtr[13] );
            // *(gPtr + 851968) = rPtr_3[13];
            *(gPtr + 851968) = rPtr[13];
        
            // rPtr_3[14] = *(gPtr + 917504);
            // turboFFT_ZADD(rPtr_3[14], rPtr_3[14], rPtr[14] );
            // *(gPtr + 917504) = rPtr_3[14];
            *(gPtr + 917504) = rPtr[14];
        
            // rPtr_3[15] = *(gPtr + 983040);
            // turboFFT_ZADD(rPtr_3[15], rPtr_3[15], rPtr[15] );
            // *(gPtr + 983040) = rPtr_3[15];
            *(gPtr + 983040) = rPtr[15];
        
            // rPtr_3[16] = *(gPtr + 1048576);
            // turboFFT_ZADD(rPtr_3[16], rPtr_3[16], rPtr[16] );
            // *(gPtr + 1048576) = rPtr_3[16];
            *(gPtr + 1048576) = rPtr[16];
        
            // rPtr_3[17] = *(gPtr + 1114112);
            // turboFFT_ZADD(rPtr_3[17], rPtr_3[17], rPtr[17] );
            // *(gPtr + 1114112) = rPtr_3[17];
            *(gPtr + 1114112) = rPtr[17];
        
            // rPtr_3[18] = *(gPtr + 1179648);
            // turboFFT_ZADD(rPtr_3[18], rPtr_3[18], rPtr[18] );
            // *(gPtr + 1179648) = rPtr_3[18];
            *(gPtr + 1179648) = rPtr[18];
        
            // rPtr_3[19] = *(gPtr + 1245184);
            // turboFFT_ZADD(rPtr_3[19], rPtr_3[19], rPtr[19] );
            // *(gPtr + 1245184) = rPtr_3[19];
            *(gPtr + 1245184) = rPtr[19];
        
            // rPtr_3[20] = *(gPtr + 1310720);
            // turboFFT_ZADD(rPtr_3[20], rPtr_3[20], rPtr[20] );
            // *(gPtr + 1310720) = rPtr_3[20];
            *(gPtr + 1310720) = rPtr[20];
        
            // rPtr_3[21] = *(gPtr + 1376256);
            // turboFFT_ZADD(rPtr_3[21], rPtr_3[21], rPtr[21] );
            // *(gPtr + 1376256) = rPtr_3[21];
            *(gPtr + 1376256) = rPtr[21];
        
            // rPtr_3[22] = *(gPtr + 1441792);
            // turboFFT_ZADD(rPtr_3[22], rPtr_3[22], rPtr[22] );
            // *(gPtr + 1441792) = rPtr_3[22];
            *(gPtr + 1441792) = rPtr[22];
        
            // rPtr_3[23] = *(gPtr + 1507328);
            // turboFFT_ZADD(rPtr_3[23], rPtr_3[23], rPtr[23] );
            // *(gPtr + 1507328) = rPtr_3[23];
            *(gPtr + 1507328) = rPtr[23];
        
            // rPtr_3[24] = *(gPtr + 1572864);
            // turboFFT_ZADD(rPtr_3[24], rPtr_3[24], rPtr[24] );
            // *(gPtr + 1572864) = rPtr_3[24];
            *(gPtr + 1572864) = rPtr[24];
        
            // rPtr_3[25] = *(gPtr + 1638400);
            // turboFFT_ZADD(rPtr_3[25], rPtr_3[25], rPtr[25] );
            // *(gPtr + 1638400) = rPtr_3[25];
            *(gPtr + 1638400) = rPtr[25];
        
            // rPtr_3[26] = *(gPtr + 1703936);
            // turboFFT_ZADD(rPtr_3[26], rPtr_3[26], rPtr[26] );
            // *(gPtr + 1703936) = rPtr_3[26];
            *(gPtr + 1703936) = rPtr[26];
        
            // rPtr_3[27] = *(gPtr + 1769472);
            // turboFFT_ZADD(rPtr_3[27], rPtr_3[27], rPtr[27] );
            // *(gPtr + 1769472) = rPtr_3[27];
            *(gPtr + 1769472) = rPtr[27];
        
            // rPtr_3[28] = *(gPtr + 1835008);
            // turboFFT_ZADD(rPtr_3[28], rPtr_3[28], rPtr[28] );
            // *(gPtr + 1835008) = rPtr_3[28];
            *(gPtr + 1835008) = rPtr[28];
        
            // rPtr_3[29] = *(gPtr + 1900544);
            // turboFFT_ZADD(rPtr_3[29], rPtr_3[29], rPtr[29] );
            // *(gPtr + 1900544) = rPtr_3[29];
            *(gPtr + 1900544) = rPtr[29];
        
            // rPtr_3[30] = *(gPtr + 1966080);
            // turboFFT_ZADD(rPtr_3[30], rPtr_3[30], rPtr[30] );
            // *(gPtr + 1966080) = rPtr_3[30];
            *(gPtr + 1966080) = rPtr[30];
        
            // rPtr_3[31] = *(gPtr + 2031616);
            // turboFFT_ZADD(rPtr_3[31], rPtr_3[31], rPtr[31] );
            // *(gPtr + 2031616) = rPtr_3[31];
            *(gPtr + 2031616) = rPtr[31];
        }}