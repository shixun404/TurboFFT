
#include "../../../TurboFFT_radix_2_template.h"
template<>
__global__ void fft_radix_2<double2, 28, 2, 0, 0>(double2* inputs, double2* outputs, double2* twiddle, double2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    double2* shared = (double2*) ext_shared;
    int threadblock_per_SM = 2;
    int tb_gap = threadblock_per_SM * 108;
    int delta_bid = ((blockIdx.x / tb_gap) ==  (gridDim.x / tb_gap)) ? (gridDim.x % tb_gap) : tb_gap;
    double2 r[3];
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
    double2* gPtr;
    double2* shPtr;
    double2 rPtr[32];
    double2 rPtr_2[32];
    double2 rPtr_3[32];
    double2 rPtr_4[32];
    double2 tmp;
    double2 tmp_1;
    double2 tmp_2;
    double2 tmp_3;
    double2 angle;
    double2 delta_angle;
    j = 0;
    k = -1;
    global_j = 0;
    global_k = 0;
    data_id = 0;
    bs_id = 0;
    shared_offset_bs = 0;
    shared_offset_data = 0;
    bx = gridDim.x - blockIdx.x - 1;
    tx = threadIdx.x;
    offset = 0;
    gPtr = inputs;
    shPtr = shared;
    
    __syncthreads();
    int bid = 0;
    for(bid = (blockIdx.x / tb_gap) * tb_gap * thread_bs + blockIdx.x % tb_gap;
                bid_cnt < thread_bs && bid < (268435456 * BS + 8192 - 1) / 8192; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 8 * 1;
    
    gPtr += (bx % 1) * 1024 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 1024;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 524288;
    bx = bx / 64;
    
    gPtr += tx % 8 * 524288;
    
    gPtr += (bx % BS * 268435456);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        rPtr[1] = *(gPtr + 32);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
        rPtr[2] = *(gPtr + 64);
        rPtr_3[2].x += rPtr[2].x;
        rPtr_3[2].y += rPtr[2].y;
        
        rPtr[3] = *(gPtr + 96);
        rPtr_3[3].x += rPtr[3].x;
        rPtr_3[3].y += rPtr[3].y;
        
        rPtr[4] = *(gPtr + 128);
        rPtr_3[4].x += rPtr[4].x;
        rPtr_3[4].y += rPtr[4].y;
        
        rPtr[5] = *(gPtr + 160);
        rPtr_3[5].x += rPtr[5].x;
        rPtr_3[5].y += rPtr[5].y;
        
        rPtr[6] = *(gPtr + 192);
        rPtr_3[6].x += rPtr[6].x;
        rPtr_3[6].y += rPtr[6].y;
        
        rPtr[7] = *(gPtr + 224);
        rPtr_3[7].x += rPtr[7].x;
        rPtr_3[7].y += rPtr[7].y;
        
        rPtr[8] = *(gPtr + 256);
        rPtr_3[8].x += rPtr[8].x;
        rPtr_3[8].y += rPtr[8].y;
        
        rPtr[9] = *(gPtr + 288);
        rPtr_3[9].x += rPtr[9].x;
        rPtr_3[9].y += rPtr[9].y;
        
        rPtr[10] = *(gPtr + 320);
        rPtr_3[10].x += rPtr[10].x;
        rPtr_3[10].y += rPtr[10].y;
        
        rPtr[11] = *(gPtr + 352);
        rPtr_3[11].x += rPtr[11].x;
        rPtr_3[11].y += rPtr[11].y;
        
        rPtr[12] = *(gPtr + 384);
        rPtr_3[12].x += rPtr[12].x;
        rPtr_3[12].y += rPtr[12].y;
        
        rPtr[13] = *(gPtr + 416);
        rPtr_3[13].x += rPtr[13].x;
        rPtr_3[13].y += rPtr[13].y;
        
        rPtr[14] = *(gPtr + 448);
        rPtr_3[14].x += rPtr[14].x;
        rPtr_3[14].y += rPtr[14].y;
        
        rPtr[15] = *(gPtr + 480);
        rPtr_3[15].x += rPtr[15].x;
        rPtr_3[15].y += rPtr[15].y;
        
        rPtr[16] = *(gPtr + 512);
        rPtr_3[16].x += rPtr[16].x;
        rPtr_3[16].y += rPtr[16].y;
        
        rPtr[17] = *(gPtr + 544);
        rPtr_3[17].x += rPtr[17].x;
        rPtr_3[17].y += rPtr[17].y;
        
        rPtr[18] = *(gPtr + 576);
        rPtr_3[18].x += rPtr[18].x;
        rPtr_3[18].y += rPtr[18].y;
        
        rPtr[19] = *(gPtr + 608);
        rPtr_3[19].x += rPtr[19].x;
        rPtr_3[19].y += rPtr[19].y;
        
        rPtr[20] = *(gPtr + 640);
        rPtr_3[20].x += rPtr[20].x;
        rPtr_3[20].y += rPtr[20].y;
        
        rPtr[21] = *(gPtr + 672);
        rPtr_3[21].x += rPtr[21].x;
        rPtr_3[21].y += rPtr[21].y;
        
        rPtr[22] = *(gPtr + 704);
        rPtr_3[22].x += rPtr[22].x;
        rPtr_3[22].y += rPtr[22].y;
        
        rPtr[23] = *(gPtr + 736);
        rPtr_3[23].x += rPtr[23].x;
        rPtr_3[23].y += rPtr[23].y;
        
        rPtr[24] = *(gPtr + 768);
        rPtr_3[24].x += rPtr[24].x;
        rPtr_3[24].y += rPtr[24].y;
        
        rPtr[25] = *(gPtr + 800);
        rPtr_3[25].x += rPtr[25].x;
        rPtr_3[25].y += rPtr[25].y;
        
        rPtr[26] = *(gPtr + 832);
        rPtr_3[26].x += rPtr[26].x;
        rPtr_3[26].y += rPtr[26].y;
        
        rPtr[27] = *(gPtr + 864);
        rPtr_3[27].x += rPtr[27].x;
        rPtr_3[27].y += rPtr[27].y;
        
        rPtr[28] = *(gPtr + 896);
        rPtr_3[28].x += rPtr[28].x;
        rPtr_3[28].y += rPtr[28].y;
        
        rPtr[29] = *(gPtr + 928);
        rPtr_3[29].x += rPtr[29].x;
        rPtr_3[29].y += rPtr[29].y;
        
        rPtr[30] = *(gPtr + 960);
        rPtr_3[30].x += rPtr[30].x;
        rPtr_3[30].y += rPtr[30].y;
        
        rPtr[31] = *(gPtr + 992);
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
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    delta_angle = twiddle[1023 + j];
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
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    gPtr += tx / 8 * 262144;
    
    gPtr += (bx % 1) * 1024 * 262144;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 512;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 1;
    bx = bx / 64;
    
    gPtr += tx % 8 * 1;
    
    gPtr += (bx % BS * 268435456);
    
            *(gPtr + 0) = rPtr[0];
            rPtr_4[0].x += rPtr[0].x;
            rPtr_4[0].y += rPtr[0].y;
            
            *(gPtr + 8388608) = rPtr[16];
            rPtr_4[1].x += rPtr[16].x;
            rPtr_4[1].y += rPtr[16].y;
            
            *(gPtr + 16777216) = rPtr[8];
            rPtr_4[2].x += rPtr[8].x;
            rPtr_4[2].y += rPtr[8].y;
            
            *(gPtr + 25165824) = rPtr[24];
            rPtr_4[3].x += rPtr[24].x;
            rPtr_4[3].y += rPtr[24].y;
            
            *(gPtr + 33554432) = rPtr[4];
            rPtr_4[4].x += rPtr[4].x;
            rPtr_4[4].y += rPtr[4].y;
            
            *(gPtr + 41943040) = rPtr[20];
            rPtr_4[5].x += rPtr[20].x;
            rPtr_4[5].y += rPtr[20].y;
            
            *(gPtr + 50331648) = rPtr[12];
            rPtr_4[6].x += rPtr[12].x;
            rPtr_4[6].y += rPtr[12].y;
            
            *(gPtr + 58720256) = rPtr[28];
            rPtr_4[7].x += rPtr[28].x;
            rPtr_4[7].y += rPtr[28].y;
            
            *(gPtr + 67108864) = rPtr[2];
            rPtr_4[8].x += rPtr[2].x;
            rPtr_4[8].y += rPtr[2].y;
            
            *(gPtr + 75497472) = rPtr[18];
            rPtr_4[9].x += rPtr[18].x;
            rPtr_4[9].y += rPtr[18].y;
            
            *(gPtr + 83886080) = rPtr[10];
            rPtr_4[10].x += rPtr[10].x;
            rPtr_4[10].y += rPtr[10].y;
            
            *(gPtr + 92274688) = rPtr[26];
            rPtr_4[11].x += rPtr[26].x;
            rPtr_4[11].y += rPtr[26].y;
            
            *(gPtr + 100663296) = rPtr[6];
            rPtr_4[12].x += rPtr[6].x;
            rPtr_4[12].y += rPtr[6].y;
            
            *(gPtr + 109051904) = rPtr[22];
            rPtr_4[13].x += rPtr[22].x;
            rPtr_4[13].y += rPtr[22].y;
            
            *(gPtr + 117440512) = rPtr[14];
            rPtr_4[14].x += rPtr[14].x;
            rPtr_4[14].y += rPtr[14].y;
            
            *(gPtr + 125829120) = rPtr[30];
            rPtr_4[15].x += rPtr[30].x;
            rPtr_4[15].y += rPtr[30].y;
            
            *(gPtr + 134217728) = rPtr[1];
            rPtr_4[16].x += rPtr[1].x;
            rPtr_4[16].y += rPtr[1].y;
            
            *(gPtr + 142606336) = rPtr[17];
            rPtr_4[17].x += rPtr[17].x;
            rPtr_4[17].y += rPtr[17].y;
            
            *(gPtr + 150994944) = rPtr[9];
            rPtr_4[18].x += rPtr[9].x;
            rPtr_4[18].y += rPtr[9].y;
            
            *(gPtr + 159383552) = rPtr[25];
            rPtr_4[19].x += rPtr[25].x;
            rPtr_4[19].y += rPtr[25].y;
            
            *(gPtr + 167772160) = rPtr[5];
            rPtr_4[20].x += rPtr[5].x;
            rPtr_4[20].y += rPtr[5].y;
            
            *(gPtr + 176160768) = rPtr[21];
            rPtr_4[21].x += rPtr[21].x;
            rPtr_4[21].y += rPtr[21].y;
            
            *(gPtr + 184549376) = rPtr[13];
            rPtr_4[22].x += rPtr[13].x;
            rPtr_4[22].y += rPtr[13].y;
            
            *(gPtr + 192937984) = rPtr[29];
            rPtr_4[23].x += rPtr[29].x;
            rPtr_4[23].y += rPtr[29].y;
            
            *(gPtr + 201326592) = rPtr[3];
            rPtr_4[24].x += rPtr[3].x;
            rPtr_4[24].y += rPtr[3].y;
            
            *(gPtr + 209715200) = rPtr[19];
            rPtr_4[25].x += rPtr[19].x;
            rPtr_4[25].y += rPtr[19].y;
            
            *(gPtr + 218103808) = rPtr[11];
            rPtr_4[26].x += rPtr[11].x;
            rPtr_4[26].y += rPtr[11].y;
            
            *(gPtr + 226492416) = rPtr[27];
            rPtr_4[27].x += rPtr[27].x;
            rPtr_4[27].y += rPtr[27].y;
            
            *(gPtr + 234881024) = rPtr[7];
            rPtr_4[28].x += rPtr[7].x;
            rPtr_4[28].y += rPtr[7].y;
            
            *(gPtr + 243269632) = rPtr[23];
            rPtr_4[29].x += rPtr[23].x;
            rPtr_4[29].y += rPtr[23].y;
            
            *(gPtr + 251658240) = rPtr[15];
            rPtr_4[30].x += rPtr[15].x;
            rPtr_4[30].y += rPtr[15].y;
            
            *(gPtr + 260046848) = rPtr[31];
            rPtr_4[31].x += rPtr[31].x;
            rPtr_4[31].y += rPtr[31].y;
            
    }
    
}

#include "../../../TurboFFT_radix_2_template.h"
template<>
__global__ void fft_radix_2<double2, 28, 2, 1, 0>(double2* inputs, double2* outputs, double2* twiddle, double2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    double2* shared = (double2*) ext_shared;
    int threadblock_per_SM = 2;
    int tb_gap = threadblock_per_SM * 108;
    int delta_bid = ((blockIdx.x / tb_gap) ==  (gridDim.x / tb_gap)) ? (gridDim.x % tb_gap) : tb_gap;
    double2 r[3];
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
    double2* gPtr;
    double2* shPtr;
    double2 rPtr[32];
    double2 rPtr_2[32];
    double2 rPtr_3[32];
    double2 rPtr_4[32];
    double2 tmp;
    double2 tmp_1;
    double2 tmp_2;
    double2 tmp_3;
    double2 angle;
    double2 delta_angle;
    j = 0;
    k = -1;
    global_j = 0;
    global_k = 0;
    data_id = 0;
    bs_id = 0;
    shared_offset_bs = 0;
    shared_offset_data = 0;
    bx = gridDim.x - blockIdx.x - 1;
    tx = threadIdx.x;
    offset = 0;
    gPtr = inputs;
    shPtr = shared;
    
    rPtr_2[0] = *(checksum_DFT + 1024 - 2 + tx + 0);
    shPtr[tx + 0] = rPtr_2[0];
    
    rPtr_2[1] = *(checksum_DFT + 1024 - 2 + tx + 256);
    shPtr[tx + 256] = rPtr_2[1];
    
    rPtr_2[2] = *(checksum_DFT + 1024 - 2 + tx + 512);
    shPtr[tx + 512] = rPtr_2[2];
    
    rPtr_2[3] = *(checksum_DFT + 1024 - 2 + tx + 768);
    shPtr[tx + 768] = rPtr_2[3];
    
    __syncthreads();
    tmp_1.x = 0;
    tmp_1.y = 0;
    tmp_2.x = 0;
    tmp_2.y = 0;
    tmp_3.x = 0;
    tmp_3.y = 0;
    
    rPtr_2[0] = *(shPtr +  tx / 8 + 0);
    rPtr_3[0].x = 0; rPtr_3[0].y = 0;
    rPtr_4[0].x = 0; rPtr_4[0].y = 0;
    
    rPtr_2[1] = *(shPtr +  tx / 8 + 32);
    rPtr_3[1].x = 0; rPtr_3[1].y = 0;
    rPtr_4[1].x = 0; rPtr_4[1].y = 0;
    
    rPtr_2[2] = *(shPtr +  tx / 8 + 64);
    rPtr_3[2].x = 0; rPtr_3[2].y = 0;
    rPtr_4[2].x = 0; rPtr_4[2].y = 0;
    
    rPtr_2[3] = *(shPtr +  tx / 8 + 96);
    rPtr_3[3].x = 0; rPtr_3[3].y = 0;
    rPtr_4[3].x = 0; rPtr_4[3].y = 0;
    
    rPtr_2[4] = *(shPtr +  tx / 8 + 128);
    rPtr_3[4].x = 0; rPtr_3[4].y = 0;
    rPtr_4[4].x = 0; rPtr_4[4].y = 0;
    
    rPtr_2[5] = *(shPtr +  tx / 8 + 160);
    rPtr_3[5].x = 0; rPtr_3[5].y = 0;
    rPtr_4[5].x = 0; rPtr_4[5].y = 0;
    
    rPtr_2[6] = *(shPtr +  tx / 8 + 192);
    rPtr_3[6].x = 0; rPtr_3[6].y = 0;
    rPtr_4[6].x = 0; rPtr_4[6].y = 0;
    
    rPtr_2[7] = *(shPtr +  tx / 8 + 224);
    rPtr_3[7].x = 0; rPtr_3[7].y = 0;
    rPtr_4[7].x = 0; rPtr_4[7].y = 0;
    
    rPtr_2[8] = *(shPtr +  tx / 8 + 256);
    rPtr_3[8].x = 0; rPtr_3[8].y = 0;
    rPtr_4[8].x = 0; rPtr_4[8].y = 0;
    
    rPtr_2[9] = *(shPtr +  tx / 8 + 288);
    rPtr_3[9].x = 0; rPtr_3[9].y = 0;
    rPtr_4[9].x = 0; rPtr_4[9].y = 0;
    
    rPtr_2[10] = *(shPtr +  tx / 8 + 320);
    rPtr_3[10].x = 0; rPtr_3[10].y = 0;
    rPtr_4[10].x = 0; rPtr_4[10].y = 0;
    
    rPtr_2[11] = *(shPtr +  tx / 8 + 352);
    rPtr_3[11].x = 0; rPtr_3[11].y = 0;
    rPtr_4[11].x = 0; rPtr_4[11].y = 0;
    
    rPtr_2[12] = *(shPtr +  tx / 8 + 384);
    rPtr_3[12].x = 0; rPtr_3[12].y = 0;
    rPtr_4[12].x = 0; rPtr_4[12].y = 0;
    
    rPtr_2[13] = *(shPtr +  tx / 8 + 416);
    rPtr_3[13].x = 0; rPtr_3[13].y = 0;
    rPtr_4[13].x = 0; rPtr_4[13].y = 0;
    
    rPtr_2[14] = *(shPtr +  tx / 8 + 448);
    rPtr_3[14].x = 0; rPtr_3[14].y = 0;
    rPtr_4[14].x = 0; rPtr_4[14].y = 0;
    
    rPtr_2[15] = *(shPtr +  tx / 8 + 480);
    rPtr_3[15].x = 0; rPtr_3[15].y = 0;
    rPtr_4[15].x = 0; rPtr_4[15].y = 0;
    
    rPtr_2[16] = *(shPtr +  tx / 8 + 512);
    rPtr_3[16].x = 0; rPtr_3[16].y = 0;
    rPtr_4[16].x = 0; rPtr_4[16].y = 0;
    
    rPtr_2[17] = *(shPtr +  tx / 8 + 544);
    rPtr_3[17].x = 0; rPtr_3[17].y = 0;
    rPtr_4[17].x = 0; rPtr_4[17].y = 0;
    
    rPtr_2[18] = *(shPtr +  tx / 8 + 576);
    rPtr_3[18].x = 0; rPtr_3[18].y = 0;
    rPtr_4[18].x = 0; rPtr_4[18].y = 0;
    
    rPtr_2[19] = *(shPtr +  tx / 8 + 608);
    rPtr_3[19].x = 0; rPtr_3[19].y = 0;
    rPtr_4[19].x = 0; rPtr_4[19].y = 0;
    
    rPtr_2[20] = *(shPtr +  tx / 8 + 640);
    rPtr_3[20].x = 0; rPtr_3[20].y = 0;
    rPtr_4[20].x = 0; rPtr_4[20].y = 0;
    
    rPtr_2[21] = *(shPtr +  tx / 8 + 672);
    rPtr_3[21].x = 0; rPtr_3[21].y = 0;
    rPtr_4[21].x = 0; rPtr_4[21].y = 0;
    
    rPtr_2[22] = *(shPtr +  tx / 8 + 704);
    rPtr_3[22].x = 0; rPtr_3[22].y = 0;
    rPtr_4[22].x = 0; rPtr_4[22].y = 0;
    
    rPtr_2[23] = *(shPtr +  tx / 8 + 736);
    rPtr_3[23].x = 0; rPtr_3[23].y = 0;
    rPtr_4[23].x = 0; rPtr_4[23].y = 0;
    
    rPtr_2[24] = *(shPtr +  tx / 8 + 768);
    rPtr_3[24].x = 0; rPtr_3[24].y = 0;
    rPtr_4[24].x = 0; rPtr_4[24].y = 0;
    
    rPtr_2[25] = *(shPtr +  tx / 8 + 800);
    rPtr_3[25].x = 0; rPtr_3[25].y = 0;
    rPtr_4[25].x = 0; rPtr_4[25].y = 0;
    
    rPtr_2[26] = *(shPtr +  tx / 8 + 832);
    rPtr_3[26].x = 0; rPtr_3[26].y = 0;
    rPtr_4[26].x = 0; rPtr_4[26].y = 0;
    
    rPtr_2[27] = *(shPtr +  tx / 8 + 864);
    rPtr_3[27].x = 0; rPtr_3[27].y = 0;
    rPtr_4[27].x = 0; rPtr_4[27].y = 0;
    
    rPtr_2[28] = *(shPtr +  tx / 8 + 896);
    rPtr_3[28].x = 0; rPtr_3[28].y = 0;
    rPtr_4[28].x = 0; rPtr_4[28].y = 0;
    
    rPtr_2[29] = *(shPtr +  tx / 8 + 928);
    rPtr_3[29].x = 0; rPtr_3[29].y = 0;
    rPtr_4[29].x = 0; rPtr_4[29].y = 0;
    
    rPtr_2[30] = *(shPtr +  tx / 8 + 960);
    rPtr_3[30].x = 0; rPtr_3[30].y = 0;
    rPtr_4[30].x = 0; rPtr_4[30].y = 0;
    
    rPtr_2[31] = *(shPtr +  tx / 8 + 992);
    rPtr_3[31].x = 0; rPtr_3[31].y = 0;
    rPtr_4[31].x = 0; rPtr_4[31].y = 0;
    
    __syncthreads();
    int bid = 0;
    for(bid = (blockIdx.x / tb_gap) * tb_gap * thread_bs + blockIdx.x % tb_gap;
                bid_cnt < thread_bs && bid < (268435456 * BS + 8192 - 1) / 8192; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 8 * 1;
    
    gPtr += (bx % 1) * 1024 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 1024;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 524288;
    bx = bx / 64;
    
    gPtr += tx % 8 * 524288;
    
    gPtr += (bx % BS * 268435456);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        // tmp = checksum_DFT[tx / 8 + 0];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[0], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[0], rPtr_2[0])
        turboFFT_ZMUL(tmp, rPtr[0], rPtr_2[0])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[1] = *(gPtr + 32);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
        // tmp = checksum_DFT[tx / 8 + 32];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[1], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[1], rPtr_2[1])
        turboFFT_ZMUL(tmp, rPtr[1], rPtr_2[1])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[2] = *(gPtr + 64);
        rPtr_3[2].x += rPtr[2].x;
        rPtr_3[2].y += rPtr[2].y;
        
        // tmp = checksum_DFT[tx / 8 + 64];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[2], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[2], rPtr_2[2])
        turboFFT_ZMUL(tmp, rPtr[2], rPtr_2[2])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[3] = *(gPtr + 96);
        rPtr_3[3].x += rPtr[3].x;
        rPtr_3[3].y += rPtr[3].y;
        
        // tmp = checksum_DFT[tx / 8 + 96];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[3], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[3], rPtr_2[3])
        turboFFT_ZMUL(tmp, rPtr[3], rPtr_2[3])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[4] = *(gPtr + 128);
        rPtr_3[4].x += rPtr[4].x;
        rPtr_3[4].y += rPtr[4].y;
        
        // tmp = checksum_DFT[tx / 8 + 128];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[4], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[4], rPtr_2[4])
        turboFFT_ZMUL(tmp, rPtr[4], rPtr_2[4])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[5] = *(gPtr + 160);
        rPtr_3[5].x += rPtr[5].x;
        rPtr_3[5].y += rPtr[5].y;
        
        // tmp = checksum_DFT[tx / 8 + 160];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[5], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[5], rPtr_2[5])
        turboFFT_ZMUL(tmp, rPtr[5], rPtr_2[5])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[6] = *(gPtr + 192);
        rPtr_3[6].x += rPtr[6].x;
        rPtr_3[6].y += rPtr[6].y;
        
        // tmp = checksum_DFT[tx / 8 + 192];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[6], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[6], rPtr_2[6])
        turboFFT_ZMUL(tmp, rPtr[6], rPtr_2[6])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[7] = *(gPtr + 224);
        rPtr_3[7].x += rPtr[7].x;
        rPtr_3[7].y += rPtr[7].y;
        
        // tmp = checksum_DFT[tx / 8 + 224];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[7], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[7], rPtr_2[7])
        turboFFT_ZMUL(tmp, rPtr[7], rPtr_2[7])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[8] = *(gPtr + 256);
        rPtr_3[8].x += rPtr[8].x;
        rPtr_3[8].y += rPtr[8].y;
        
        // tmp = checksum_DFT[tx / 8 + 256];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[8], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[8], rPtr_2[8])
        turboFFT_ZMUL(tmp, rPtr[8], rPtr_2[8])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[9] = *(gPtr + 288);
        rPtr_3[9].x += rPtr[9].x;
        rPtr_3[9].y += rPtr[9].y;
        
        // tmp = checksum_DFT[tx / 8 + 288];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[9], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[9], rPtr_2[9])
        turboFFT_ZMUL(tmp, rPtr[9], rPtr_2[9])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[10] = *(gPtr + 320);
        rPtr_3[10].x += rPtr[10].x;
        rPtr_3[10].y += rPtr[10].y;
        
        // tmp = checksum_DFT[tx / 8 + 320];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[10], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[10], rPtr_2[10])
        turboFFT_ZMUL(tmp, rPtr[10], rPtr_2[10])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[11] = *(gPtr + 352);
        rPtr_3[11].x += rPtr[11].x;
        rPtr_3[11].y += rPtr[11].y;
        
        // tmp = checksum_DFT[tx / 8 + 352];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[11], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[11], rPtr_2[11])
        turboFFT_ZMUL(tmp, rPtr[11], rPtr_2[11])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[12] = *(gPtr + 384);
        rPtr_3[12].x += rPtr[12].x;
        rPtr_3[12].y += rPtr[12].y;
        
        // tmp = checksum_DFT[tx / 8 + 384];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[12], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[12], rPtr_2[12])
        turboFFT_ZMUL(tmp, rPtr[12], rPtr_2[12])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[13] = *(gPtr + 416);
        rPtr_3[13].x += rPtr[13].x;
        rPtr_3[13].y += rPtr[13].y;
        
        // tmp = checksum_DFT[tx / 8 + 416];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[13], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[13], rPtr_2[13])
        turboFFT_ZMUL(tmp, rPtr[13], rPtr_2[13])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[14] = *(gPtr + 448);
        rPtr_3[14].x += rPtr[14].x;
        rPtr_3[14].y += rPtr[14].y;
        
        // tmp = checksum_DFT[tx / 8 + 448];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[14], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[14], rPtr_2[14])
        turboFFT_ZMUL(tmp, rPtr[14], rPtr_2[14])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[15] = *(gPtr + 480);
        rPtr_3[15].x += rPtr[15].x;
        rPtr_3[15].y += rPtr[15].y;
        
        // tmp = checksum_DFT[tx / 8 + 480];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[15], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[15], rPtr_2[15])
        turboFFT_ZMUL(tmp, rPtr[15], rPtr_2[15])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[16] = *(gPtr + 512);
        rPtr_3[16].x += rPtr[16].x;
        rPtr_3[16].y += rPtr[16].y;
        
        // tmp = checksum_DFT[tx / 8 + 512];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[16], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[16], rPtr_2[16])
        turboFFT_ZMUL(tmp, rPtr[16], rPtr_2[16])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[17] = *(gPtr + 544);
        rPtr_3[17].x += rPtr[17].x;
        rPtr_3[17].y += rPtr[17].y;
        
        // tmp = checksum_DFT[tx / 8 + 544];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[17], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[17], rPtr_2[17])
        turboFFT_ZMUL(tmp, rPtr[17], rPtr_2[17])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[18] = *(gPtr + 576);
        rPtr_3[18].x += rPtr[18].x;
        rPtr_3[18].y += rPtr[18].y;
        
        // tmp = checksum_DFT[tx / 8 + 576];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[18], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[18], rPtr_2[18])
        turboFFT_ZMUL(tmp, rPtr[18], rPtr_2[18])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[19] = *(gPtr + 608);
        rPtr_3[19].x += rPtr[19].x;
        rPtr_3[19].y += rPtr[19].y;
        
        // tmp = checksum_DFT[tx / 8 + 608];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[19], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[19], rPtr_2[19])
        turboFFT_ZMUL(tmp, rPtr[19], rPtr_2[19])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[20] = *(gPtr + 640);
        rPtr_3[20].x += rPtr[20].x;
        rPtr_3[20].y += rPtr[20].y;
        
        // tmp = checksum_DFT[tx / 8 + 640];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[20], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[20], rPtr_2[20])
        turboFFT_ZMUL(tmp, rPtr[20], rPtr_2[20])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[21] = *(gPtr + 672);
        rPtr_3[21].x += rPtr[21].x;
        rPtr_3[21].y += rPtr[21].y;
        
        // tmp = checksum_DFT[tx / 8 + 672];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[21], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[21], rPtr_2[21])
        turboFFT_ZMUL(tmp, rPtr[21], rPtr_2[21])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[22] = *(gPtr + 704);
        rPtr_3[22].x += rPtr[22].x;
        rPtr_3[22].y += rPtr[22].y;
        
        // tmp = checksum_DFT[tx / 8 + 704];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[22], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[22], rPtr_2[22])
        turboFFT_ZMUL(tmp, rPtr[22], rPtr_2[22])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[23] = *(gPtr + 736);
        rPtr_3[23].x += rPtr[23].x;
        rPtr_3[23].y += rPtr[23].y;
        
        // tmp = checksum_DFT[tx / 8 + 736];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[23], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[23], rPtr_2[23])
        turboFFT_ZMUL(tmp, rPtr[23], rPtr_2[23])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[24] = *(gPtr + 768);
        rPtr_3[24].x += rPtr[24].x;
        rPtr_3[24].y += rPtr[24].y;
        
        // tmp = checksum_DFT[tx / 8 + 768];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[24], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[24], rPtr_2[24])
        turboFFT_ZMUL(tmp, rPtr[24], rPtr_2[24])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[25] = *(gPtr + 800);
        rPtr_3[25].x += rPtr[25].x;
        rPtr_3[25].y += rPtr[25].y;
        
        // tmp = checksum_DFT[tx / 8 + 800];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[25], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[25], rPtr_2[25])
        turboFFT_ZMUL(tmp, rPtr[25], rPtr_2[25])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[26] = *(gPtr + 832);
        rPtr_3[26].x += rPtr[26].x;
        rPtr_3[26].y += rPtr[26].y;
        
        // tmp = checksum_DFT[tx / 8 + 832];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[26], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[26], rPtr_2[26])
        turboFFT_ZMUL(tmp, rPtr[26], rPtr_2[26])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[27] = *(gPtr + 864);
        rPtr_3[27].x += rPtr[27].x;
        rPtr_3[27].y += rPtr[27].y;
        
        // tmp = checksum_DFT[tx / 8 + 864];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[27], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[27], rPtr_2[27])
        turboFFT_ZMUL(tmp, rPtr[27], rPtr_2[27])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[28] = *(gPtr + 896);
        rPtr_3[28].x += rPtr[28].x;
        rPtr_3[28].y += rPtr[28].y;
        
        // tmp = checksum_DFT[tx / 8 + 896];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[28], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[28], rPtr_2[28])
        turboFFT_ZMUL(tmp, rPtr[28], rPtr_2[28])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[29] = *(gPtr + 928);
        rPtr_3[29].x += rPtr[29].x;
        rPtr_3[29].y += rPtr[29].y;
        
        // tmp = checksum_DFT[tx / 8 + 928];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[29], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[29], rPtr_2[29])
        turboFFT_ZMUL(tmp, rPtr[29], rPtr_2[29])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[30] = *(gPtr + 960);
        rPtr_3[30].x += rPtr[30].x;
        rPtr_3[30].y += rPtr[30].y;
        
        // tmp = checksum_DFT[tx / 8 + 960];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[30], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[30], rPtr_2[30])
        turboFFT_ZMUL(tmp, rPtr[30], rPtr_2[30])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[31] = *(gPtr + 992);
        rPtr_3[31].x += rPtr[31].x;
        rPtr_3[31].y += rPtr[31].y;
        
        // tmp = checksum_DFT[tx / 8 + 992];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[31], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[31], rPtr_2[31])
        turboFFT_ZMUL(tmp, rPtr[31], rPtr_2[31])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        // tmp_3.x += bid_cnt * (rPtr[0].x + rPtr[0].y) * 1024;
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[16]);
    turboFFT_ZSUB(rPtr[16], tmp, rPtr[16]);
    tmp = rPtr[16];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[17]);
    turboFFT_ZSUB(rPtr[17], tmp, rPtr[17]);
    tmp = rPtr[17];
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    delta_angle = twiddle[1023 + j];
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
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    gPtr += tx / 8 * 262144;
    
    gPtr += (bx % 1) * 1024 * 262144;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 512;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 1;
    bx = bx / 64;
    
    gPtr += tx % 8 * 1;
    
    gPtr += (bx % BS * 268435456);
    
        // 1's vector
        // tmp_3.y -=  (rPtr[0].y + rPtr[0].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[0].y + rPtr[0].x);
        turboFFT_ZMUL(tmp, rPtr[0],r[(0 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[0], r[(0 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[0], r[(0 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[16].y + rPtr[16].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[16].y + rPtr[16].x);
        turboFFT_ZMUL(tmp, rPtr[16],r[(32 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[16], r[(32 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[16], r[(32 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[8].y + rPtr[8].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[8].y + rPtr[8].x);
        turboFFT_ZMUL(tmp, rPtr[8],r[(64 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[8], r[(64 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[8], r[(64 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[24].y + rPtr[24].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[24].y + rPtr[24].x);
        turboFFT_ZMUL(tmp, rPtr[24],r[(96 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[24], r[(96 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[24], r[(96 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[4].y + rPtr[4].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[4].y + rPtr[4].x);
        turboFFT_ZMUL(tmp, rPtr[4],r[(128 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[4], r[(128 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[4], r[(128 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[20].y + rPtr[20].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[20].y + rPtr[20].x);
        turboFFT_ZMUL(tmp, rPtr[20],r[(160 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[20], r[(160 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[20], r[(160 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[12].y + rPtr[12].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[12].y + rPtr[12].x);
        turboFFT_ZMUL(tmp, rPtr[12],r[(192 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[12], r[(192 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[12], r[(192 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[28].y + rPtr[28].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[28].y + rPtr[28].x);
        turboFFT_ZMUL(tmp, rPtr[28],r[(224 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[28], r[(224 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[28], r[(224 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[2].y + rPtr[2].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[2].y + rPtr[2].x);
        turboFFT_ZMUL(tmp, rPtr[2],r[(256 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[2], r[(256 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[2], r[(256 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[18].y + rPtr[18].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[18].y + rPtr[18].x);
        turboFFT_ZMUL(tmp, rPtr[18],r[(288 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[18], r[(288 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[18], r[(288 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[10].y + rPtr[10].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[10].y + rPtr[10].x);
        turboFFT_ZMUL(tmp, rPtr[10],r[(320 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[10], r[(320 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[10], r[(320 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[26].y + rPtr[26].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[26].y + rPtr[26].x);
        turboFFT_ZMUL(tmp, rPtr[26],r[(352 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[26], r[(352 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[26], r[(352 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[6].y + rPtr[6].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[6].y + rPtr[6].x);
        turboFFT_ZMUL(tmp, rPtr[6],r[(384 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[6], r[(384 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[6], r[(384 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[22].y + rPtr[22].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[22].y + rPtr[22].x);
        turboFFT_ZMUL(tmp, rPtr[22],r[(416 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[22], r[(416 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[22], r[(416 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[14].y + rPtr[14].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[14].y + rPtr[14].x);
        turboFFT_ZMUL(tmp, rPtr[14],r[(448 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[14], r[(448 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[14], r[(448 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[30].y + rPtr[30].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[30].y + rPtr[30].x);
        turboFFT_ZMUL(tmp, rPtr[30],r[(480 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[30], r[(480 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[30], r[(480 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[1].y + rPtr[1].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[1].y + rPtr[1].x);
        turboFFT_ZMUL(tmp, rPtr[1],r[(512 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[1], r[(512 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[1], r[(512 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[17].y + rPtr[17].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[17].y + rPtr[17].x);
        turboFFT_ZMUL(tmp, rPtr[17],r[(544 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[17], r[(544 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[17], r[(544 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[9].y + rPtr[9].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[9].y + rPtr[9].x);
        turboFFT_ZMUL(tmp, rPtr[9],r[(576 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[9], r[(576 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[9], r[(576 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[25].y + rPtr[25].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[25].y + rPtr[25].x);
        turboFFT_ZMUL(tmp, rPtr[25],r[(608 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[25], r[(608 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[25], r[(608 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[5].y + rPtr[5].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[5].y + rPtr[5].x);
        turboFFT_ZMUL(tmp, rPtr[5],r[(640 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[5], r[(640 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[5], r[(640 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[21].y + rPtr[21].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[21].y + rPtr[21].x);
        turboFFT_ZMUL(tmp, rPtr[21],r[(672 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[21], r[(672 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[21], r[(672 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[13].y + rPtr[13].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[13].y + rPtr[13].x);
        turboFFT_ZMUL(tmp, rPtr[13],r[(704 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[13], r[(704 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[13], r[(704 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[29].y + rPtr[29].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[29].y + rPtr[29].x);
        turboFFT_ZMUL(tmp, rPtr[29],r[(736 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[29], r[(736 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[29], r[(736 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[3].y + rPtr[3].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[3].y + rPtr[3].x);
        turboFFT_ZMUL(tmp, rPtr[3],r[(768 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[3], r[(768 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[3], r[(768 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[19].y + rPtr[19].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[19].y + rPtr[19].x);
        turboFFT_ZMUL(tmp, rPtr[19],r[(800 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[19], r[(800 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[19], r[(800 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[11].y + rPtr[11].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[11].y + rPtr[11].x);
        turboFFT_ZMUL(tmp, rPtr[11],r[(832 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[11], r[(832 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[11], r[(832 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[27].y + rPtr[27].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[27].y + rPtr[27].x);
        turboFFT_ZMUL(tmp, rPtr[27],r[(864 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[27], r[(864 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[27], r[(864 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[7].y + rPtr[7].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[7].y + rPtr[7].x);
        turboFFT_ZMUL(tmp, rPtr[7],r[(896 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[7], r[(896 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[7], r[(896 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[23].y + rPtr[23].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[23].y + rPtr[23].x);
        turboFFT_ZMUL(tmp, rPtr[23],r[(928 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[23], r[(928 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[23], r[(928 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[15].y + rPtr[15].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[15].y + rPtr[15].x);
        turboFFT_ZMUL(tmp, rPtr[15],r[(960 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[15], r[(960 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[15], r[(960 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[31].y + rPtr[31].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[31].y + rPtr[31].x);
        turboFFT_ZMUL(tmp, rPtr[31],r[(992 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[31], r[(992 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[31], r[(992 + tx / 8) % 3])
        
            *(gPtr + 0) = rPtr[0];
            rPtr_4[0].x += rPtr[0].x;
            rPtr_4[0].y += rPtr[0].y;
            
            *(gPtr + 8388608) = rPtr[16];
            rPtr_4[1].x += rPtr[16].x;
            rPtr_4[1].y += rPtr[16].y;
            
            *(gPtr + 16777216) = rPtr[8];
            rPtr_4[2].x += rPtr[8].x;
            rPtr_4[2].y += rPtr[8].y;
            
            *(gPtr + 25165824) = rPtr[24];
            rPtr_4[3].x += rPtr[24].x;
            rPtr_4[3].y += rPtr[24].y;
            
            *(gPtr + 33554432) = rPtr[4];
            rPtr_4[4].x += rPtr[4].x;
            rPtr_4[4].y += rPtr[4].y;
            
            *(gPtr + 41943040) = rPtr[20];
            rPtr_4[5].x += rPtr[20].x;
            rPtr_4[5].y += rPtr[20].y;
            
            *(gPtr + 50331648) = rPtr[12];
            rPtr_4[6].x += rPtr[12].x;
            rPtr_4[6].y += rPtr[12].y;
            
            *(gPtr + 58720256) = rPtr[28];
            rPtr_4[7].x += rPtr[28].x;
            rPtr_4[7].y += rPtr[28].y;
            
            *(gPtr + 67108864) = rPtr[2];
            rPtr_4[8].x += rPtr[2].x;
            rPtr_4[8].y += rPtr[2].y;
            
            *(gPtr + 75497472) = rPtr[18];
            rPtr_4[9].x += rPtr[18].x;
            rPtr_4[9].y += rPtr[18].y;
            
            *(gPtr + 83886080) = rPtr[10];
            rPtr_4[10].x += rPtr[10].x;
            rPtr_4[10].y += rPtr[10].y;
            
            *(gPtr + 92274688) = rPtr[26];
            rPtr_4[11].x += rPtr[26].x;
            rPtr_4[11].y += rPtr[26].y;
            
            *(gPtr + 100663296) = rPtr[6];
            rPtr_4[12].x += rPtr[6].x;
            rPtr_4[12].y += rPtr[6].y;
            
            *(gPtr + 109051904) = rPtr[22];
            rPtr_4[13].x += rPtr[22].x;
            rPtr_4[13].y += rPtr[22].y;
            
            *(gPtr + 117440512) = rPtr[14];
            rPtr_4[14].x += rPtr[14].x;
            rPtr_4[14].y += rPtr[14].y;
            
            *(gPtr + 125829120) = rPtr[30];
            rPtr_4[15].x += rPtr[30].x;
            rPtr_4[15].y += rPtr[30].y;
            
            *(gPtr + 134217728) = rPtr[1];
            rPtr_4[16].x += rPtr[1].x;
            rPtr_4[16].y += rPtr[1].y;
            
            *(gPtr + 142606336) = rPtr[17];
            rPtr_4[17].x += rPtr[17].x;
            rPtr_4[17].y += rPtr[17].y;
            
            *(gPtr + 150994944) = rPtr[9];
            rPtr_4[18].x += rPtr[9].x;
            rPtr_4[18].y += rPtr[9].y;
            
            *(gPtr + 159383552) = rPtr[25];
            rPtr_4[19].x += rPtr[25].x;
            rPtr_4[19].y += rPtr[25].y;
            
            *(gPtr + 167772160) = rPtr[5];
            rPtr_4[20].x += rPtr[5].x;
            rPtr_4[20].y += rPtr[5].y;
            
            *(gPtr + 176160768) = rPtr[21];
            rPtr_4[21].x += rPtr[21].x;
            rPtr_4[21].y += rPtr[21].y;
            
            *(gPtr + 184549376) = rPtr[13];
            rPtr_4[22].x += rPtr[13].x;
            rPtr_4[22].y += rPtr[13].y;
            
            *(gPtr + 192937984) = rPtr[29];
            rPtr_4[23].x += rPtr[29].x;
            rPtr_4[23].y += rPtr[29].y;
            
            *(gPtr + 201326592) = rPtr[3];
            rPtr_4[24].x += rPtr[3].x;
            rPtr_4[24].y += rPtr[3].y;
            
            *(gPtr + 209715200) = rPtr[19];
            rPtr_4[25].x += rPtr[19].x;
            rPtr_4[25].y += rPtr[19].y;
            
            *(gPtr + 218103808) = rPtr[11];
            rPtr_4[26].x += rPtr[11].x;
            rPtr_4[26].y += rPtr[11].y;
            
            *(gPtr + 226492416) = rPtr[27];
            rPtr_4[27].x += rPtr[27].x;
            rPtr_4[27].y += rPtr[27].y;
            
            *(gPtr + 234881024) = rPtr[7];
            rPtr_4[28].x += rPtr[7].x;
            rPtr_4[28].y += rPtr[7].y;
            
            *(gPtr + 243269632) = rPtr[23];
            rPtr_4[29].x += rPtr[23].x;
            rPtr_4[29].y += rPtr[23].y;
            
            *(gPtr + 251658240) = rPtr[15];
            rPtr_4[30].x += rPtr[15].x;
            rPtr_4[30].y += rPtr[15].y;
            
            *(gPtr + 260046848) = rPtr[31];
            rPtr_4[31].x += rPtr[31].x;
            rPtr_4[31].y += rPtr[31].y;
            
        if(bid_cnt==thread_bs)
        
        {
        
        // 1's vector
        // tmp.x = (tx / 8 == 0) ? (rPtr_3[0].y + rPtr_3[0].x) * 1024: 0;
        // tmp.y = (tx / 8 == 0) ? (abs(rPtr_3[0].y) + abs(rPtr_3[0].x)) * 1024: 0;
        tmp = tmp_1;
        tmp_1.y += tmp.x;
        tmp_1.x = (abs(tmp.y) + abs(tmp.x));
        
        // 1's vector
        // tmp.x = (tx / 8 == 0) ? tmp_3.x : 0;
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
        
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(tx == 0)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if((blockIdx.x % thread_bs + 1) != round(abs(tmp_3.y) / abs(tmp_1.y)) && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 )  printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            //                                         bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x, tmp_3.y, tmp_3.y / tmp_1.y);
            // if(abs(tmp_1.y / tmp_1.x) > 0.001)printf("2, bid=%d bx=%d, by=%d, tx=%d: %f, %f, %f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
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
__global__ void fft_radix_2<double2, 28, 2, 1, 1>(double2* inputs, double2* outputs, double2* twiddle, double2* checksum_DFT, int BS, int thread_bs) {
    int bid_cnt = 0;
    
    double2* shared = (double2*) ext_shared;
    int threadblock_per_SM = 2;
    int tb_gap = threadblock_per_SM * 108;
    int delta_bid = ((blockIdx.x / tb_gap) ==  (gridDim.x / tb_gap)) ? (gridDim.x % tb_gap) : tb_gap;
    double2 r[3];
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
    double2* gPtr;
    double2* shPtr;
    double2 rPtr[32];
    double2 rPtr_2[32];
    double2 rPtr_3[32];
    double2 rPtr_4[32];
    double2 tmp;
    double2 tmp_1;
    double2 tmp_2;
    double2 tmp_3;
    double2 angle;
    double2 delta_angle;
    j = 0;
    k = -1;
    global_j = 0;
    global_k = 0;
    data_id = 0;
    bs_id = 0;
    shared_offset_bs = 0;
    shared_offset_data = 0;
    bx = gridDim.x - blockIdx.x - 1;
    tx = threadIdx.x;
    offset = 0;
    gPtr = inputs;
    shPtr = shared;
    
    rPtr_2[0] = *(checksum_DFT + 1024 - 2 + tx + 0);
    shPtr[tx + 0] = rPtr_2[0];
    
    rPtr_2[1] = *(checksum_DFT + 1024 - 2 + tx + 256);
    shPtr[tx + 256] = rPtr_2[1];
    
    rPtr_2[2] = *(checksum_DFT + 1024 - 2 + tx + 512);
    shPtr[tx + 512] = rPtr_2[2];
    
    rPtr_2[3] = *(checksum_DFT + 1024 - 2 + tx + 768);
    shPtr[tx + 768] = rPtr_2[3];
    
    __syncthreads();
    tmp_1.x = 0;
    tmp_1.y = 0;
    tmp_2.x = 0;
    tmp_2.y = 0;
    tmp_3.x = 0;
    tmp_3.y = 0;
    
    rPtr_2[0] = *(shPtr +  tx / 8 + 0);
    rPtr_3[0].x = 0; rPtr_3[0].y = 0;
    rPtr_4[0].x = 0; rPtr_4[0].y = 0;
    
    rPtr_2[1] = *(shPtr +  tx / 8 + 32);
    rPtr_3[1].x = 0; rPtr_3[1].y = 0;
    rPtr_4[1].x = 0; rPtr_4[1].y = 0;
    
    rPtr_2[2] = *(shPtr +  tx / 8 + 64);
    rPtr_3[2].x = 0; rPtr_3[2].y = 0;
    rPtr_4[2].x = 0; rPtr_4[2].y = 0;
    
    rPtr_2[3] = *(shPtr +  tx / 8 + 96);
    rPtr_3[3].x = 0; rPtr_3[3].y = 0;
    rPtr_4[3].x = 0; rPtr_4[3].y = 0;
    
    rPtr_2[4] = *(shPtr +  tx / 8 + 128);
    rPtr_3[4].x = 0; rPtr_3[4].y = 0;
    rPtr_4[4].x = 0; rPtr_4[4].y = 0;
    
    rPtr_2[5] = *(shPtr +  tx / 8 + 160);
    rPtr_3[5].x = 0; rPtr_3[5].y = 0;
    rPtr_4[5].x = 0; rPtr_4[5].y = 0;
    
    rPtr_2[6] = *(shPtr +  tx / 8 + 192);
    rPtr_3[6].x = 0; rPtr_3[6].y = 0;
    rPtr_4[6].x = 0; rPtr_4[6].y = 0;
    
    rPtr_2[7] = *(shPtr +  tx / 8 + 224);
    rPtr_3[7].x = 0; rPtr_3[7].y = 0;
    rPtr_4[7].x = 0; rPtr_4[7].y = 0;
    
    rPtr_2[8] = *(shPtr +  tx / 8 + 256);
    rPtr_3[8].x = 0; rPtr_3[8].y = 0;
    rPtr_4[8].x = 0; rPtr_4[8].y = 0;
    
    rPtr_2[9] = *(shPtr +  tx / 8 + 288);
    rPtr_3[9].x = 0; rPtr_3[9].y = 0;
    rPtr_4[9].x = 0; rPtr_4[9].y = 0;
    
    rPtr_2[10] = *(shPtr +  tx / 8 + 320);
    rPtr_3[10].x = 0; rPtr_3[10].y = 0;
    rPtr_4[10].x = 0; rPtr_4[10].y = 0;
    
    rPtr_2[11] = *(shPtr +  tx / 8 + 352);
    rPtr_3[11].x = 0; rPtr_3[11].y = 0;
    rPtr_4[11].x = 0; rPtr_4[11].y = 0;
    
    rPtr_2[12] = *(shPtr +  tx / 8 + 384);
    rPtr_3[12].x = 0; rPtr_3[12].y = 0;
    rPtr_4[12].x = 0; rPtr_4[12].y = 0;
    
    rPtr_2[13] = *(shPtr +  tx / 8 + 416);
    rPtr_3[13].x = 0; rPtr_3[13].y = 0;
    rPtr_4[13].x = 0; rPtr_4[13].y = 0;
    
    rPtr_2[14] = *(shPtr +  tx / 8 + 448);
    rPtr_3[14].x = 0; rPtr_3[14].y = 0;
    rPtr_4[14].x = 0; rPtr_4[14].y = 0;
    
    rPtr_2[15] = *(shPtr +  tx / 8 + 480);
    rPtr_3[15].x = 0; rPtr_3[15].y = 0;
    rPtr_4[15].x = 0; rPtr_4[15].y = 0;
    
    rPtr_2[16] = *(shPtr +  tx / 8 + 512);
    rPtr_3[16].x = 0; rPtr_3[16].y = 0;
    rPtr_4[16].x = 0; rPtr_4[16].y = 0;
    
    rPtr_2[17] = *(shPtr +  tx / 8 + 544);
    rPtr_3[17].x = 0; rPtr_3[17].y = 0;
    rPtr_4[17].x = 0; rPtr_4[17].y = 0;
    
    rPtr_2[18] = *(shPtr +  tx / 8 + 576);
    rPtr_3[18].x = 0; rPtr_3[18].y = 0;
    rPtr_4[18].x = 0; rPtr_4[18].y = 0;
    
    rPtr_2[19] = *(shPtr +  tx / 8 + 608);
    rPtr_3[19].x = 0; rPtr_3[19].y = 0;
    rPtr_4[19].x = 0; rPtr_4[19].y = 0;
    
    rPtr_2[20] = *(shPtr +  tx / 8 + 640);
    rPtr_3[20].x = 0; rPtr_3[20].y = 0;
    rPtr_4[20].x = 0; rPtr_4[20].y = 0;
    
    rPtr_2[21] = *(shPtr +  tx / 8 + 672);
    rPtr_3[21].x = 0; rPtr_3[21].y = 0;
    rPtr_4[21].x = 0; rPtr_4[21].y = 0;
    
    rPtr_2[22] = *(shPtr +  tx / 8 + 704);
    rPtr_3[22].x = 0; rPtr_3[22].y = 0;
    rPtr_4[22].x = 0; rPtr_4[22].y = 0;
    
    rPtr_2[23] = *(shPtr +  tx / 8 + 736);
    rPtr_3[23].x = 0; rPtr_3[23].y = 0;
    rPtr_4[23].x = 0; rPtr_4[23].y = 0;
    
    rPtr_2[24] = *(shPtr +  tx / 8 + 768);
    rPtr_3[24].x = 0; rPtr_3[24].y = 0;
    rPtr_4[24].x = 0; rPtr_4[24].y = 0;
    
    rPtr_2[25] = *(shPtr +  tx / 8 + 800);
    rPtr_3[25].x = 0; rPtr_3[25].y = 0;
    rPtr_4[25].x = 0; rPtr_4[25].y = 0;
    
    rPtr_2[26] = *(shPtr +  tx / 8 + 832);
    rPtr_3[26].x = 0; rPtr_3[26].y = 0;
    rPtr_4[26].x = 0; rPtr_4[26].y = 0;
    
    rPtr_2[27] = *(shPtr +  tx / 8 + 864);
    rPtr_3[27].x = 0; rPtr_3[27].y = 0;
    rPtr_4[27].x = 0; rPtr_4[27].y = 0;
    
    rPtr_2[28] = *(shPtr +  tx / 8 + 896);
    rPtr_3[28].x = 0; rPtr_3[28].y = 0;
    rPtr_4[28].x = 0; rPtr_4[28].y = 0;
    
    rPtr_2[29] = *(shPtr +  tx / 8 + 928);
    rPtr_3[29].x = 0; rPtr_3[29].y = 0;
    rPtr_4[29].x = 0; rPtr_4[29].y = 0;
    
    rPtr_2[30] = *(shPtr +  tx / 8 + 960);
    rPtr_3[30].x = 0; rPtr_3[30].y = 0;
    rPtr_4[30].x = 0; rPtr_4[30].y = 0;
    
    rPtr_2[31] = *(shPtr +  tx / 8 + 992);
    rPtr_3[31].x = 0; rPtr_3[31].y = 0;
    rPtr_4[31].x = 0; rPtr_4[31].y = 0;
    
    __syncthreads();
    int bid = 0;
    for(bid = (blockIdx.x / tb_gap) * tb_gap * thread_bs + blockIdx.x % tb_gap;
                bid_cnt < thread_bs && bid < (268435456 * BS + 8192 - 1) / 8192; bid += delta_bid)
    {
    bid_cnt += 1;
            
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 8 * 1;
    
    gPtr += (bx % 1) * 1024 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 1024;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 524288;
    bx = bx / 64;
    
    gPtr += tx % 8 * 524288;
    
    gPtr += (bx % BS * 268435456);
    
        rPtr[0] = *(gPtr + 0);
        rPtr_3[0].x += rPtr[0].x;
        rPtr_3[0].y += rPtr[0].y;
        
        // tmp = checksum_DFT[tx / 8 + 0];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[0], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[0], rPtr_2[0])
        turboFFT_ZMUL(tmp, rPtr[0], rPtr_2[0])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[1] = *(gPtr + 32);
        rPtr_3[1].x += rPtr[1].x;
        rPtr_3[1].y += rPtr[1].y;
        
        // tmp = checksum_DFT[tx / 8 + 32];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[1], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[1], rPtr_2[1])
        turboFFT_ZMUL(tmp, rPtr[1], rPtr_2[1])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[2] = *(gPtr + 64);
        rPtr_3[2].x += rPtr[2].x;
        rPtr_3[2].y += rPtr[2].y;
        
        // tmp = checksum_DFT[tx / 8 + 64];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[2], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[2], rPtr_2[2])
        turboFFT_ZMUL(tmp, rPtr[2], rPtr_2[2])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[3] = *(gPtr + 96);
        rPtr_3[3].x += rPtr[3].x;
        rPtr_3[3].y += rPtr[3].y;
        
        // tmp = checksum_DFT[tx / 8 + 96];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[3], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[3], rPtr_2[3])
        turboFFT_ZMUL(tmp, rPtr[3], rPtr_2[3])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[4] = *(gPtr + 128);
        rPtr_3[4].x += rPtr[4].x;
        rPtr_3[4].y += rPtr[4].y;
        
        // tmp = checksum_DFT[tx / 8 + 128];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[4], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[4], rPtr_2[4])
        turboFFT_ZMUL(tmp, rPtr[4], rPtr_2[4])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[5] = *(gPtr + 160);
        rPtr_3[5].x += rPtr[5].x;
        rPtr_3[5].y += rPtr[5].y;
        
        // tmp = checksum_DFT[tx / 8 + 160];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[5], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[5], rPtr_2[5])
        turboFFT_ZMUL(tmp, rPtr[5], rPtr_2[5])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[6] = *(gPtr + 192);
        rPtr_3[6].x += rPtr[6].x;
        rPtr_3[6].y += rPtr[6].y;
        
        // tmp = checksum_DFT[tx / 8 + 192];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[6], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[6], rPtr_2[6])
        turboFFT_ZMUL(tmp, rPtr[6], rPtr_2[6])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[7] = *(gPtr + 224);
        rPtr_3[7].x += rPtr[7].x;
        rPtr_3[7].y += rPtr[7].y;
        
        // tmp = checksum_DFT[tx / 8 + 224];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[7], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[7], rPtr_2[7])
        turboFFT_ZMUL(tmp, rPtr[7], rPtr_2[7])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[8] = *(gPtr + 256);
        rPtr_3[8].x += rPtr[8].x;
        rPtr_3[8].y += rPtr[8].y;
        
        // tmp = checksum_DFT[tx / 8 + 256];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[8], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[8], rPtr_2[8])
        turboFFT_ZMUL(tmp, rPtr[8], rPtr_2[8])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[9] = *(gPtr + 288);
        rPtr_3[9].x += rPtr[9].x;
        rPtr_3[9].y += rPtr[9].y;
        
        // tmp = checksum_DFT[tx / 8 + 288];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[9], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[9], rPtr_2[9])
        turboFFT_ZMUL(tmp, rPtr[9], rPtr_2[9])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[10] = *(gPtr + 320);
        rPtr_3[10].x += rPtr[10].x;
        rPtr_3[10].y += rPtr[10].y;
        
        // tmp = checksum_DFT[tx / 8 + 320];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[10], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[10], rPtr_2[10])
        turboFFT_ZMUL(tmp, rPtr[10], rPtr_2[10])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[11] = *(gPtr + 352);
        rPtr_3[11].x += rPtr[11].x;
        rPtr_3[11].y += rPtr[11].y;
        
        // tmp = checksum_DFT[tx / 8 + 352];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[11], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[11], rPtr_2[11])
        turboFFT_ZMUL(tmp, rPtr[11], rPtr_2[11])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[12] = *(gPtr + 384);
        rPtr_3[12].x += rPtr[12].x;
        rPtr_3[12].y += rPtr[12].y;
        
        // tmp = checksum_DFT[tx / 8 + 384];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[12], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[12], rPtr_2[12])
        turboFFT_ZMUL(tmp, rPtr[12], rPtr_2[12])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[13] = *(gPtr + 416);
        rPtr_3[13].x += rPtr[13].x;
        rPtr_3[13].y += rPtr[13].y;
        
        // tmp = checksum_DFT[tx / 8 + 416];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[13], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[13], rPtr_2[13])
        turboFFT_ZMUL(tmp, rPtr[13], rPtr_2[13])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[14] = *(gPtr + 448);
        rPtr_3[14].x += rPtr[14].x;
        rPtr_3[14].y += rPtr[14].y;
        
        // tmp = checksum_DFT[tx / 8 + 448];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[14], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[14], rPtr_2[14])
        turboFFT_ZMUL(tmp, rPtr[14], rPtr_2[14])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[15] = *(gPtr + 480);
        rPtr_3[15].x += rPtr[15].x;
        rPtr_3[15].y += rPtr[15].y;
        
        // tmp = checksum_DFT[tx / 8 + 480];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[15], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[15], rPtr_2[15])
        turboFFT_ZMUL(tmp, rPtr[15], rPtr_2[15])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[16] = *(gPtr + 512);
        rPtr_3[16].x += rPtr[16].x;
        rPtr_3[16].y += rPtr[16].y;
        
        // tmp = checksum_DFT[tx / 8 + 512];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[16], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[16], rPtr_2[16])
        turboFFT_ZMUL(tmp, rPtr[16], rPtr_2[16])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[17] = *(gPtr + 544);
        rPtr_3[17].x += rPtr[17].x;
        rPtr_3[17].y += rPtr[17].y;
        
        // tmp = checksum_DFT[tx / 8 + 544];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[17], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[17], rPtr_2[17])
        turboFFT_ZMUL(tmp, rPtr[17], rPtr_2[17])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[18] = *(gPtr + 576);
        rPtr_3[18].x += rPtr[18].x;
        rPtr_3[18].y += rPtr[18].y;
        
        // tmp = checksum_DFT[tx / 8 + 576];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[18], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[18], rPtr_2[18])
        turboFFT_ZMUL(tmp, rPtr[18], rPtr_2[18])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[19] = *(gPtr + 608);
        rPtr_3[19].x += rPtr[19].x;
        rPtr_3[19].y += rPtr[19].y;
        
        // tmp = checksum_DFT[tx / 8 + 608];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[19], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[19], rPtr_2[19])
        turboFFT_ZMUL(tmp, rPtr[19], rPtr_2[19])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[20] = *(gPtr + 640);
        rPtr_3[20].x += rPtr[20].x;
        rPtr_3[20].y += rPtr[20].y;
        
        // tmp = checksum_DFT[tx / 8 + 640];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[20], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[20], rPtr_2[20])
        turboFFT_ZMUL(tmp, rPtr[20], rPtr_2[20])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[21] = *(gPtr + 672);
        rPtr_3[21].x += rPtr[21].x;
        rPtr_3[21].y += rPtr[21].y;
        
        // tmp = checksum_DFT[tx / 8 + 672];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[21], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[21], rPtr_2[21])
        turboFFT_ZMUL(tmp, rPtr[21], rPtr_2[21])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[22] = *(gPtr + 704);
        rPtr_3[22].x += rPtr[22].x;
        rPtr_3[22].y += rPtr[22].y;
        
        // tmp = checksum_DFT[tx / 8 + 704];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[22], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[22], rPtr_2[22])
        turboFFT_ZMUL(tmp, rPtr[22], rPtr_2[22])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[23] = *(gPtr + 736);
        rPtr_3[23].x += rPtr[23].x;
        rPtr_3[23].y += rPtr[23].y;
        
        // tmp = checksum_DFT[tx / 8 + 736];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[23], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[23], rPtr_2[23])
        turboFFT_ZMUL(tmp, rPtr[23], rPtr_2[23])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[24] = *(gPtr + 768);
        rPtr_3[24].x += rPtr[24].x;
        rPtr_3[24].y += rPtr[24].y;
        
        // tmp = checksum_DFT[tx / 8 + 768];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[24], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[24], rPtr_2[24])
        turboFFT_ZMUL(tmp, rPtr[24], rPtr_2[24])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[25] = *(gPtr + 800);
        rPtr_3[25].x += rPtr[25].x;
        rPtr_3[25].y += rPtr[25].y;
        
        // tmp = checksum_DFT[tx / 8 + 800];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[25], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[25], rPtr_2[25])
        turboFFT_ZMUL(tmp, rPtr[25], rPtr_2[25])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[26] = *(gPtr + 832);
        rPtr_3[26].x += rPtr[26].x;
        rPtr_3[26].y += rPtr[26].y;
        
        // tmp = checksum_DFT[tx / 8 + 832];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[26], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[26], rPtr_2[26])
        turboFFT_ZMUL(tmp, rPtr[26], rPtr_2[26])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[27] = *(gPtr + 864);
        rPtr_3[27].x += rPtr[27].x;
        rPtr_3[27].y += rPtr[27].y;
        
        // tmp = checksum_DFT[tx / 8 + 864];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[27], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[27], rPtr_2[27])
        turboFFT_ZMUL(tmp, rPtr[27], rPtr_2[27])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[28] = *(gPtr + 896);
        rPtr_3[28].x += rPtr[28].x;
        rPtr_3[28].y += rPtr[28].y;
        
        // tmp = checksum_DFT[tx / 8 + 896];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[28], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[28], rPtr_2[28])
        turboFFT_ZMUL(tmp, rPtr[28], rPtr_2[28])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[29] = *(gPtr + 928);
        rPtr_3[29].x += rPtr[29].x;
        rPtr_3[29].y += rPtr[29].y;
        
        // tmp = checksum_DFT[tx / 8 + 928];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[29], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[29], rPtr_2[29])
        turboFFT_ZMUL(tmp, rPtr[29], rPtr_2[29])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[30] = *(gPtr + 960);
        rPtr_3[30].x += rPtr[30].x;
        rPtr_3[30].y += rPtr[30].y;
        
        // tmp = checksum_DFT[tx / 8 + 960];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[30], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[30], rPtr_2[30])
        turboFFT_ZMUL(tmp, rPtr[30], rPtr_2[30])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        rPtr[31] = *(gPtr + 992);
        rPtr_3[31].x += rPtr[31].x;
        rPtr_3[31].y += rPtr[31].y;
        
        // tmp = checksum_DFT[tx / 8 + 992];
        // turboFFT_ZMUL_ACC(tmp_1, rPtr[31], tmp);
        //  turboFFT_ZMUL_ACC(tmp_1, rPtr[31], rPtr_2[31])
        turboFFT_ZMUL(tmp, rPtr[31], rPtr_2[31])
        tmp_1.x += (tmp.x + tmp.y);
        tmp_3.x += bid_cnt * (tmp.x + tmp.y);
        
        // tmp_3.x += bid_cnt * (rPtr[0].x + rPtr[0].y) * 1024;
        
        rPtr[0].x += (threadIdx.x == 0 && bid_cnt == (blockIdx.x % thread_bs + 1)) ? 100: 0;
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[16]);
    turboFFT_ZSUB(rPtr[16], tmp, rPtr[16]);
    tmp = rPtr[16];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[17]);
    turboFFT_ZSUB(rPtr[17], tmp, rPtr[17]);
    tmp = rPtr[17];
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    delta_angle = twiddle[1023 + j];
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
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    gPtr += tx / 8 * 262144;
    
    gPtr += (bx % 1) * 1024 * 262144;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 512;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 1;
    bx = bx / 64;
    
    gPtr += tx % 8 * 1;
    
    gPtr += (bx % BS * 268435456);
    
        // 1's vector
        // tmp_3.y -=  (rPtr[0].y + rPtr[0].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[0].y + rPtr[0].x);
        turboFFT_ZMUL(tmp, rPtr[0],r[(0 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[0], r[(0 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[0], r[(0 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[16].y + rPtr[16].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[16].y + rPtr[16].x);
        turboFFT_ZMUL(tmp, rPtr[16],r[(32 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[16], r[(32 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[16], r[(32 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[8].y + rPtr[8].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[8].y + rPtr[8].x);
        turboFFT_ZMUL(tmp, rPtr[8],r[(64 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[8], r[(64 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[8], r[(64 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[24].y + rPtr[24].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[24].y + rPtr[24].x);
        turboFFT_ZMUL(tmp, rPtr[24],r[(96 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[24], r[(96 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[24], r[(96 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[4].y + rPtr[4].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[4].y + rPtr[4].x);
        turboFFT_ZMUL(tmp, rPtr[4],r[(128 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[4], r[(128 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[4], r[(128 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[20].y + rPtr[20].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[20].y + rPtr[20].x);
        turboFFT_ZMUL(tmp, rPtr[20],r[(160 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[20], r[(160 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[20], r[(160 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[12].y + rPtr[12].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[12].y + rPtr[12].x);
        turboFFT_ZMUL(tmp, rPtr[12],r[(192 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[12], r[(192 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[12], r[(192 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[28].y + rPtr[28].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[28].y + rPtr[28].x);
        turboFFT_ZMUL(tmp, rPtr[28],r[(224 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[28], r[(224 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[28], r[(224 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[2].y + rPtr[2].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[2].y + rPtr[2].x);
        turboFFT_ZMUL(tmp, rPtr[2],r[(256 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[2], r[(256 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[2], r[(256 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[18].y + rPtr[18].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[18].y + rPtr[18].x);
        turboFFT_ZMUL(tmp, rPtr[18],r[(288 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[18], r[(288 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[18], r[(288 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[10].y + rPtr[10].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[10].y + rPtr[10].x);
        turboFFT_ZMUL(tmp, rPtr[10],r[(320 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[10], r[(320 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[10], r[(320 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[26].y + rPtr[26].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[26].y + rPtr[26].x);
        turboFFT_ZMUL(tmp, rPtr[26],r[(352 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[26], r[(352 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[26], r[(352 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[6].y + rPtr[6].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[6].y + rPtr[6].x);
        turboFFT_ZMUL(tmp, rPtr[6],r[(384 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[6], r[(384 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[6], r[(384 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[22].y + rPtr[22].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[22].y + rPtr[22].x);
        turboFFT_ZMUL(tmp, rPtr[22],r[(416 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[22], r[(416 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[22], r[(416 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[14].y + rPtr[14].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[14].y + rPtr[14].x);
        turboFFT_ZMUL(tmp, rPtr[14],r[(448 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[14], r[(448 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[14], r[(448 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[30].y + rPtr[30].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[30].y + rPtr[30].x);
        turboFFT_ZMUL(tmp, rPtr[30],r[(480 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[30], r[(480 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[30], r[(480 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[1].y + rPtr[1].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[1].y + rPtr[1].x);
        turboFFT_ZMUL(tmp, rPtr[1],r[(512 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[1], r[(512 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[1], r[(512 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[17].y + rPtr[17].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[17].y + rPtr[17].x);
        turboFFT_ZMUL(tmp, rPtr[17],r[(544 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[17], r[(544 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[17], r[(544 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[9].y + rPtr[9].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[9].y + rPtr[9].x);
        turboFFT_ZMUL(tmp, rPtr[9],r[(576 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[9], r[(576 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[9], r[(576 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[25].y + rPtr[25].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[25].y + rPtr[25].x);
        turboFFT_ZMUL(tmp, rPtr[25],r[(608 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[25], r[(608 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[25], r[(608 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[5].y + rPtr[5].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[5].y + rPtr[5].x);
        turboFFT_ZMUL(tmp, rPtr[5],r[(640 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[5], r[(640 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[5], r[(640 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[21].y + rPtr[21].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[21].y + rPtr[21].x);
        turboFFT_ZMUL(tmp, rPtr[21],r[(672 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[21], r[(672 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[21], r[(672 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[13].y + rPtr[13].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[13].y + rPtr[13].x);
        turboFFT_ZMUL(tmp, rPtr[13],r[(704 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[13], r[(704 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[13], r[(704 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[29].y + rPtr[29].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[29].y + rPtr[29].x);
        turboFFT_ZMUL(tmp, rPtr[29],r[(736 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[29], r[(736 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[29], r[(736 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[3].y + rPtr[3].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[3].y + rPtr[3].x);
        turboFFT_ZMUL(tmp, rPtr[3],r[(768 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[3], r[(768 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[3], r[(768 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[19].y + rPtr[19].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[19].y + rPtr[19].x);
        turboFFT_ZMUL(tmp, rPtr[19],r[(800 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[19], r[(800 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[19], r[(800 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[11].y + rPtr[11].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[11].y + rPtr[11].x);
        turboFFT_ZMUL(tmp, rPtr[11],r[(832 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[11], r[(832 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[11], r[(832 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[27].y + rPtr[27].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[27].y + rPtr[27].x);
        turboFFT_ZMUL(tmp, rPtr[27],r[(864 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[27], r[(864 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[27], r[(864 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[7].y + rPtr[7].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[7].y + rPtr[7].x);
        turboFFT_ZMUL(tmp, rPtr[7],r[(896 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[7], r[(896 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[7], r[(896 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[23].y + rPtr[23].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[23].y + rPtr[23].x);
        turboFFT_ZMUL(tmp, rPtr[23],r[(928 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[23], r[(928 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[23], r[(928 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[15].y + rPtr[15].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[15].y + rPtr[15].x);
        turboFFT_ZMUL(tmp, rPtr[15],r[(960 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[15], r[(960 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[15], r[(960 + tx / 8) % 3])
        
        // 1's vector
        // tmp_3.y -=  (rPtr[31].y + rPtr[31].x) * bid_cnt;
        // tmp_1.y -=  (rPtr[31].y + rPtr[31].x);
        turboFFT_ZMUL(tmp, rPtr[31],r[(992 + tx / 8) % 3])
        tmp_1.y -= (tmp.x + tmp.y);
        tmp_3.y -= (tmp.y + tmp.x) * bid_cnt;
        // turboFFT_ZMUL_NACC(tmp_1,  rPtr[31], r[(992 + tx / 8) % 3])
        // turboFFT_ZMUL_NACC(tmp_3,  rPtr[31], r[(992 + tx / 8) % 3])
        
            *(gPtr + 0) = rPtr[0];
            rPtr_4[0].x += rPtr[0].x;
            rPtr_4[0].y += rPtr[0].y;
            
            *(gPtr + 8388608) = rPtr[16];
            rPtr_4[1].x += rPtr[16].x;
            rPtr_4[1].y += rPtr[16].y;
            
            *(gPtr + 16777216) = rPtr[8];
            rPtr_4[2].x += rPtr[8].x;
            rPtr_4[2].y += rPtr[8].y;
            
            *(gPtr + 25165824) = rPtr[24];
            rPtr_4[3].x += rPtr[24].x;
            rPtr_4[3].y += rPtr[24].y;
            
            *(gPtr + 33554432) = rPtr[4];
            rPtr_4[4].x += rPtr[4].x;
            rPtr_4[4].y += rPtr[4].y;
            
            *(gPtr + 41943040) = rPtr[20];
            rPtr_4[5].x += rPtr[20].x;
            rPtr_4[5].y += rPtr[20].y;
            
            *(gPtr + 50331648) = rPtr[12];
            rPtr_4[6].x += rPtr[12].x;
            rPtr_4[6].y += rPtr[12].y;
            
            *(gPtr + 58720256) = rPtr[28];
            rPtr_4[7].x += rPtr[28].x;
            rPtr_4[7].y += rPtr[28].y;
            
            *(gPtr + 67108864) = rPtr[2];
            rPtr_4[8].x += rPtr[2].x;
            rPtr_4[8].y += rPtr[2].y;
            
            *(gPtr + 75497472) = rPtr[18];
            rPtr_4[9].x += rPtr[18].x;
            rPtr_4[9].y += rPtr[18].y;
            
            *(gPtr + 83886080) = rPtr[10];
            rPtr_4[10].x += rPtr[10].x;
            rPtr_4[10].y += rPtr[10].y;
            
            *(gPtr + 92274688) = rPtr[26];
            rPtr_4[11].x += rPtr[26].x;
            rPtr_4[11].y += rPtr[26].y;
            
            *(gPtr + 100663296) = rPtr[6];
            rPtr_4[12].x += rPtr[6].x;
            rPtr_4[12].y += rPtr[6].y;
            
            *(gPtr + 109051904) = rPtr[22];
            rPtr_4[13].x += rPtr[22].x;
            rPtr_4[13].y += rPtr[22].y;
            
            *(gPtr + 117440512) = rPtr[14];
            rPtr_4[14].x += rPtr[14].x;
            rPtr_4[14].y += rPtr[14].y;
            
            *(gPtr + 125829120) = rPtr[30];
            rPtr_4[15].x += rPtr[30].x;
            rPtr_4[15].y += rPtr[30].y;
            
            *(gPtr + 134217728) = rPtr[1];
            rPtr_4[16].x += rPtr[1].x;
            rPtr_4[16].y += rPtr[1].y;
            
            *(gPtr + 142606336) = rPtr[17];
            rPtr_4[17].x += rPtr[17].x;
            rPtr_4[17].y += rPtr[17].y;
            
            *(gPtr + 150994944) = rPtr[9];
            rPtr_4[18].x += rPtr[9].x;
            rPtr_4[18].y += rPtr[9].y;
            
            *(gPtr + 159383552) = rPtr[25];
            rPtr_4[19].x += rPtr[25].x;
            rPtr_4[19].y += rPtr[25].y;
            
            *(gPtr + 167772160) = rPtr[5];
            rPtr_4[20].x += rPtr[5].x;
            rPtr_4[20].y += rPtr[5].y;
            
            *(gPtr + 176160768) = rPtr[21];
            rPtr_4[21].x += rPtr[21].x;
            rPtr_4[21].y += rPtr[21].y;
            
            *(gPtr + 184549376) = rPtr[13];
            rPtr_4[22].x += rPtr[13].x;
            rPtr_4[22].y += rPtr[13].y;
            
            *(gPtr + 192937984) = rPtr[29];
            rPtr_4[23].x += rPtr[29].x;
            rPtr_4[23].y += rPtr[29].y;
            
            *(gPtr + 201326592) = rPtr[3];
            rPtr_4[24].x += rPtr[3].x;
            rPtr_4[24].y += rPtr[3].y;
            
            *(gPtr + 209715200) = rPtr[19];
            rPtr_4[25].x += rPtr[19].x;
            rPtr_4[25].y += rPtr[19].y;
            
            *(gPtr + 218103808) = rPtr[11];
            rPtr_4[26].x += rPtr[11].x;
            rPtr_4[26].y += rPtr[11].y;
            
            *(gPtr + 226492416) = rPtr[27];
            rPtr_4[27].x += rPtr[27].x;
            rPtr_4[27].y += rPtr[27].y;
            
            *(gPtr + 234881024) = rPtr[7];
            rPtr_4[28].x += rPtr[7].x;
            rPtr_4[28].y += rPtr[7].y;
            
            *(gPtr + 243269632) = rPtr[23];
            rPtr_4[29].x += rPtr[23].x;
            rPtr_4[29].y += rPtr[23].y;
            
            *(gPtr + 251658240) = rPtr[15];
            rPtr_4[30].x += rPtr[15].x;
            rPtr_4[30].y += rPtr[15].y;
            
            *(gPtr + 260046848) = rPtr[31];
            rPtr_4[31].x += rPtr[31].x;
            rPtr_4[31].y += rPtr[31].y;
            
        if(bid_cnt==thread_bs)
        
        {
        
        // 1's vector
        // tmp.x = (tx / 8 == 0) ? (rPtr_3[0].y + rPtr_3[0].x) * 1024: 0;
        // tmp.y = (tx / 8 == 0) ? (abs(rPtr_3[0].y) + abs(rPtr_3[0].x)) * 1024: 0;
        tmp = tmp_1;
        tmp_1.y += tmp.x;
        tmp_1.x = (abs(tmp.y) + abs(tmp.x));
        
        // 1's vector
        // tmp.x = (tx / 8 == 0) ? tmp_3.x : 0;
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
        
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(tx == 0)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if(tx == 0 && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 1e-3)printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            // if((blockIdx.x % thread_bs + 1) != round(abs(tmp_3.y) / abs(tmp_1.y)) && abs(tmp_1.y) / (1000 + abs(tmp_1.x)) > 0.001 )  printf("2, bid=%d bx=%d, by=%d, tx=%d: checksum=%f, delta=%f, rel=%f, delta_3=%f, delta_3/delta=%f\n",
            //                                         bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x, tmp_3.y, tmp_3.y / tmp_1.y);
            // if(abs(tmp_1.y / tmp_1.x) > 0.001)printf("2, bid=%d bx=%d, by=%d, tx=%d: %f, %f, %f\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
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
                // if(threadIdx.x == 0)printf("bid=%d, upload=%d, bx=%d, tx=%d, k = %d\n", bid, 2, blockIdx.x, threadIdx.x, k);
                // bid = k;
                        
    bx = bid;
    tx = threadIdx.x;
    
            gPtr = inputs;
    
    gPtr += tx / 8 * 1;
    
    gPtr += (bx % 1) * 1024 * 1;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 1024;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 524288;
    bx = bx / 64;
    
    gPtr += tx % 8 * 524288;
    
    gPtr += (bx % BS * 268435456);
    
        // rPtr[0] = rPtr_3[0];
        rPtr[0] = *(gPtr + 0);
        
        // rPtr[1] = rPtr_3[1];
        rPtr[1] = *(gPtr + 32);
        
        // rPtr[2] = rPtr_3[2];
        rPtr[2] = *(gPtr + 64);
        
        // rPtr[3] = rPtr_3[3];
        rPtr[3] = *(gPtr + 96);
        
        // rPtr[4] = rPtr_3[4];
        rPtr[4] = *(gPtr + 128);
        
        // rPtr[5] = rPtr_3[5];
        rPtr[5] = *(gPtr + 160);
        
        // rPtr[6] = rPtr_3[6];
        rPtr[6] = *(gPtr + 192);
        
        // rPtr[7] = rPtr_3[7];
        rPtr[7] = *(gPtr + 224);
        
        // rPtr[8] = rPtr_3[8];
        rPtr[8] = *(gPtr + 256);
        
        // rPtr[9] = rPtr_3[9];
        rPtr[9] = *(gPtr + 288);
        
        // rPtr[10] = rPtr_3[10];
        rPtr[10] = *(gPtr + 320);
        
        // rPtr[11] = rPtr_3[11];
        rPtr[11] = *(gPtr + 352);
        
        // rPtr[12] = rPtr_3[12];
        rPtr[12] = *(gPtr + 384);
        
        // rPtr[13] = rPtr_3[13];
        rPtr[13] = *(gPtr + 416);
        
        // rPtr[14] = rPtr_3[14];
        rPtr[14] = *(gPtr + 448);
        
        // rPtr[15] = rPtr_3[15];
        rPtr[15] = *(gPtr + 480);
        
        // rPtr[16] = rPtr_3[16];
        rPtr[16] = *(gPtr + 512);
        
        // rPtr[17] = rPtr_3[17];
        rPtr[17] = *(gPtr + 544);
        
        // rPtr[18] = rPtr_3[18];
        rPtr[18] = *(gPtr + 576);
        
        // rPtr[19] = rPtr_3[19];
        rPtr[19] = *(gPtr + 608);
        
        // rPtr[20] = rPtr_3[20];
        rPtr[20] = *(gPtr + 640);
        
        // rPtr[21] = rPtr_3[21];
        rPtr[21] = *(gPtr + 672);
        
        // rPtr[22] = rPtr_3[22];
        rPtr[22] = *(gPtr + 704);
        
        // rPtr[23] = rPtr_3[23];
        rPtr[23] = *(gPtr + 736);
        
        // rPtr[24] = rPtr_3[24];
        rPtr[24] = *(gPtr + 768);
        
        // rPtr[25] = rPtr_3[25];
        rPtr[25] = *(gPtr + 800);
        
        // rPtr[26] = rPtr_3[26];
        rPtr[26] = *(gPtr + 832);
        
        // rPtr[27] = rPtr_3[27];
        rPtr[27] = *(gPtr + 864);
        
        // rPtr[28] = rPtr_3[28];
        rPtr[28] = *(gPtr + 896);
        
        // rPtr[29] = rPtr_3[29];
        rPtr[29] = *(gPtr + 928);
        
        // rPtr[30] = rPtr_3[30];
        rPtr[30] = *(gPtr + 960);
        
        // rPtr[31] = rPtr_3[31];
        rPtr[31] = *(gPtr + 992);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[16]);
    turboFFT_ZSUB(rPtr[16], tmp, rPtr[16]);
    tmp = rPtr[16];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[17]);
    turboFFT_ZSUB(rPtr[17], tmp, rPtr[17]);
    tmp = rPtr[17];
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    delta_angle = twiddle[1023 + j];
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
    
        angle.x = 0.9807852804032304;
        angle.y = -0.19509032201612825;
        turboFFT_ZMUL(rPtr[17], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[18]);
    turboFFT_ZSUB(rPtr[18], tmp, rPtr[18]);
    tmp = rPtr[18];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[18], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[19]);
    turboFFT_ZSUB(rPtr[19], tmp, rPtr[19]);
    tmp = rPtr[19];
    
        angle.x = 0.8314696123025452;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[19], tmp, angle);
        
    tmp = rPtr[4];
    turboFFT_ZADD(rPtr[4], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[20], tmp, angle);
        
    tmp = rPtr[5];
    turboFFT_ZADD(rPtr[5], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.5555702330196023;
        angle.y = -0.8314696123025452;
        turboFFT_ZMUL(rPtr[21], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[22]);
    turboFFT_ZSUB(rPtr[22], tmp, rPtr[22]);
    tmp = rPtr[22];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[22], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[23]);
    turboFFT_ZSUB(rPtr[23], tmp, rPtr[23]);
    tmp = rPtr[23];
    
        angle.x = 0.19509032201612833;
        angle.y = -0.9807852804032304;
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
    
        angle.x = -0.1950903220161282;
        angle.y = -0.9807852804032304;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[10];
    turboFFT_ZADD(rPtr[10], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[11];
    turboFFT_ZADD(rPtr[11], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = -0.555570233019602;
        angle.y = -0.8314696123025455;
        turboFFT_ZMUL(rPtr[27], tmp, angle);
        
    tmp = rPtr[12];
    turboFFT_ZADD(rPtr[12], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[28], tmp, angle);
        
    tmp = rPtr[13];
    turboFFT_ZADD(rPtr[13], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = -0.8314696123025453;
        angle.y = -0.5555702330196022;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[14];
    turboFFT_ZADD(rPtr[14], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[15];
    turboFFT_ZADD(rPtr[15], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9807852804032304;
        angle.y = -0.1950903220161286;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[8]);
    turboFFT_ZSUB(rPtr[8], tmp, rPtr[8]);
    tmp = rPtr[8];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[9]);
    turboFFT_ZSUB(rPtr[9], tmp, rPtr[9]);
    tmp = rPtr[9];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[9], tmp, angle);
        
    tmp = rPtr[2];
    turboFFT_ZADD(rPtr[2], tmp, rPtr[10]);
    turboFFT_ZSUB(rPtr[10], tmp, rPtr[10]);
    tmp = rPtr[10];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[10], tmp, angle);
        
    tmp = rPtr[3];
    turboFFT_ZADD(rPtr[3], tmp, rPtr[11]);
    turboFFT_ZSUB(rPtr[11], tmp, rPtr[11]);
    tmp = rPtr[11];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[13], tmp, angle);
        
    tmp = rPtr[6];
    turboFFT_ZADD(rPtr[6], tmp, rPtr[14]);
    turboFFT_ZSUB(rPtr[14], tmp, rPtr[14]);
    tmp = rPtr[14];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[14], tmp, angle);
        
    tmp = rPtr[7];
    turboFFT_ZADD(rPtr[7], tmp, rPtr[15]);
    turboFFT_ZSUB(rPtr[15], tmp, rPtr[15]);
    tmp = rPtr[15];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[24]);
    turboFFT_ZSUB(rPtr[24], tmp, rPtr[24]);
    tmp = rPtr[24];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[25]);
    turboFFT_ZSUB(rPtr[25], tmp, rPtr[25]);
    tmp = rPtr[25];
    
        angle.x = 0.9238795325112867;
        angle.y = -0.3826834323650898;
        turboFFT_ZMUL(rPtr[25], tmp, angle);
        
    tmp = rPtr[18];
    turboFFT_ZADD(rPtr[18], tmp, rPtr[26]);
    turboFFT_ZSUB(rPtr[26], tmp, rPtr[26]);
    tmp = rPtr[26];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
        turboFFT_ZMUL(rPtr[26], tmp, angle);
        
    tmp = rPtr[19];
    turboFFT_ZADD(rPtr[19], tmp, rPtr[27]);
    turboFFT_ZSUB(rPtr[27], tmp, rPtr[27]);
    tmp = rPtr[27];
    
        angle.x = 0.38268343236508984;
        angle.y = -0.9238795325112867;
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
    
        angle.x = -0.3826834323650897;
        angle.y = -0.9238795325112867;
        turboFFT_ZMUL(rPtr[29], tmp, angle);
        
    tmp = rPtr[22];
    turboFFT_ZADD(rPtr[22], tmp, rPtr[30]);
    turboFFT_ZSUB(rPtr[30], tmp, rPtr[30]);
    tmp = rPtr[30];
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[30], tmp, angle);
        
    tmp = rPtr[23];
    turboFFT_ZADD(rPtr[23], tmp, rPtr[31]);
    turboFFT_ZSUB(rPtr[31], tmp, rPtr[31]);
    tmp = rPtr[31];
    
        angle.x = -0.9238795325112867;
        angle.y = -0.3826834323650899;
        turboFFT_ZMUL(rPtr[31], tmp, angle);
        
    tmp = rPtr[0];
    turboFFT_ZADD(rPtr[0], tmp, rPtr[4]);
    turboFFT_ZSUB(rPtr[4], tmp, rPtr[4]);
    tmp = rPtr[4];
    
    tmp = rPtr[1];
    turboFFT_ZADD(rPtr[1], tmp, rPtr[5]);
    turboFFT_ZSUB(rPtr[5], tmp, rPtr[5]);
    tmp = rPtr[5];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[7], tmp, angle);
        
    tmp = rPtr[8];
    turboFFT_ZADD(rPtr[8], tmp, rPtr[12]);
    turboFFT_ZSUB(rPtr[12], tmp, rPtr[12]);
    tmp = rPtr[12];
    
    tmp = rPtr[9];
    turboFFT_ZADD(rPtr[9], tmp, rPtr[13]);
    turboFFT_ZSUB(rPtr[13], tmp, rPtr[13]);
    tmp = rPtr[13];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[15], tmp, angle);
        
    tmp = rPtr[16];
    turboFFT_ZADD(rPtr[16], tmp, rPtr[20]);
    turboFFT_ZSUB(rPtr[20], tmp, rPtr[20]);
    tmp = rPtr[20];
    
    tmp = rPtr[17];
    turboFFT_ZADD(rPtr[17], tmp, rPtr[21]);
    turboFFT_ZSUB(rPtr[21], tmp, rPtr[21]);
    tmp = rPtr[21];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
        turboFFT_ZMUL(rPtr[23], tmp, angle);
        
    tmp = rPtr[24];
    turboFFT_ZADD(rPtr[24], tmp, rPtr[28]);
    turboFFT_ZSUB(rPtr[28], tmp, rPtr[28]);
    tmp = rPtr[28];
    
    tmp = rPtr[25];
    turboFFT_ZADD(rPtr[25], tmp, rPtr[29]);
    turboFFT_ZSUB(rPtr[29], tmp, rPtr[29]);
    tmp = rPtr[29];
    
        angle.x = 0.7071067811865476;
        angle.y = -0.7071067811865475;
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
    
        angle.x = -0.7071067811865475;
        angle.y = -0.7071067811865476;
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
    
    gPtr += tx / 8 * 262144;
    
    gPtr += (bx % 1) * 1024 * 262144;
    bx = bx / 1;
    
    gPtr += (bx % 512) * 1 * 512;
    bx = bx / 512;
    
    gPtr += (bx % 64) * 8 * 1;
    bx = bx / 64;
    
    gPtr += tx % 8 * 1;
    
    gPtr += (bx % BS * 268435456);
    
            // turboFFT_ZSUB(rPtr[0], rPtr[0], rPtr_4[0]);
            
            // turboFFT_ZSUB(rPtr[16], rPtr[16], rPtr_4[1]);
            
            // turboFFT_ZSUB(rPtr[8], rPtr[8], rPtr_4[2]);
            
            // turboFFT_ZSUB(rPtr[24], rPtr[24], rPtr_4[3]);
            
            // turboFFT_ZSUB(rPtr[4], rPtr[4], rPtr_4[4]);
            
            // turboFFT_ZSUB(rPtr[20], rPtr[20], rPtr_4[5]);
            
            // turboFFT_ZSUB(rPtr[12], rPtr[12], rPtr_4[6]);
            
            // turboFFT_ZSUB(rPtr[28], rPtr[28], rPtr_4[7]);
            
            // turboFFT_ZSUB(rPtr[2], rPtr[2], rPtr_4[8]);
            
            // turboFFT_ZSUB(rPtr[18], rPtr[18], rPtr_4[9]);
            
            // turboFFT_ZSUB(rPtr[10], rPtr[10], rPtr_4[10]);
            
            // turboFFT_ZSUB(rPtr[26], rPtr[26], rPtr_4[11]);
            
            // turboFFT_ZSUB(rPtr[6], rPtr[6], rPtr_4[12]);
            
            // turboFFT_ZSUB(rPtr[22], rPtr[22], rPtr_4[13]);
            
            // turboFFT_ZSUB(rPtr[14], rPtr[14], rPtr_4[14]);
            
            // turboFFT_ZSUB(rPtr[30], rPtr[30], rPtr_4[15]);
            
            // turboFFT_ZSUB(rPtr[1], rPtr[1], rPtr_4[16]);
            
            // turboFFT_ZSUB(rPtr[17], rPtr[17], rPtr_4[17]);
            
            // turboFFT_ZSUB(rPtr[9], rPtr[9], rPtr_4[18]);
            
            // turboFFT_ZSUB(rPtr[25], rPtr[25], rPtr_4[19]);
            
            // turboFFT_ZSUB(rPtr[5], rPtr[5], rPtr_4[20]);
            
            // turboFFT_ZSUB(rPtr[21], rPtr[21], rPtr_4[21]);
            
            // turboFFT_ZSUB(rPtr[13], rPtr[13], rPtr_4[22]);
            
            // turboFFT_ZSUB(rPtr[29], rPtr[29], rPtr_4[23]);
            
            // turboFFT_ZSUB(rPtr[3], rPtr[3], rPtr_4[24]);
            
            // turboFFT_ZSUB(rPtr[19], rPtr[19], rPtr_4[25]);
            
            // turboFFT_ZSUB(rPtr[11], rPtr[11], rPtr_4[26]);
            
            // turboFFT_ZSUB(rPtr[27], rPtr[27], rPtr_4[27]);
            
            // turboFFT_ZSUB(rPtr[7], rPtr[7], rPtr_4[28]);
            
            // turboFFT_ZSUB(rPtr[23], rPtr[23], rPtr_4[29]);
            
            // turboFFT_ZSUB(rPtr[15], rPtr[15], rPtr_4[30]);
            
            // turboFFT_ZSUB(rPtr[31], rPtr[31], rPtr_4[31]);
            
            // rPtr_3[0] = *(gPtr + 0);
            // turboFFT_ZADD(rPtr_3[0], rPtr_3[0], rPtr[0] );
            // *(gPtr + 0) = rPtr_3[0];
            *(gPtr + 0) = rPtr[0];
        
            // rPtr_3[1] = *(gPtr + 8388608);
            // turboFFT_ZADD(rPtr_3[1], rPtr_3[1], rPtr[16] );
            // *(gPtr + 8388608) = rPtr_3[1];
            *(gPtr + 8388608) = rPtr[16];
        
            // rPtr_3[2] = *(gPtr + 16777216);
            // turboFFT_ZADD(rPtr_3[2], rPtr_3[2], rPtr[8] );
            // *(gPtr + 16777216) = rPtr_3[2];
            *(gPtr + 16777216) = rPtr[8];
        
            // rPtr_3[3] = *(gPtr + 25165824);
            // turboFFT_ZADD(rPtr_3[3], rPtr_3[3], rPtr[24] );
            // *(gPtr + 25165824) = rPtr_3[3];
            *(gPtr + 25165824) = rPtr[24];
        
            // rPtr_3[4] = *(gPtr + 33554432);
            // turboFFT_ZADD(rPtr_3[4], rPtr_3[4], rPtr[4] );
            // *(gPtr + 33554432) = rPtr_3[4];
            *(gPtr + 33554432) = rPtr[4];
        
            // rPtr_3[5] = *(gPtr + 41943040);
            // turboFFT_ZADD(rPtr_3[5], rPtr_3[5], rPtr[20] );
            // *(gPtr + 41943040) = rPtr_3[5];
            *(gPtr + 41943040) = rPtr[20];
        
            // rPtr_3[6] = *(gPtr + 50331648);
            // turboFFT_ZADD(rPtr_3[6], rPtr_3[6], rPtr[12] );
            // *(gPtr + 50331648) = rPtr_3[6];
            *(gPtr + 50331648) = rPtr[12];
        
            // rPtr_3[7] = *(gPtr + 58720256);
            // turboFFT_ZADD(rPtr_3[7], rPtr_3[7], rPtr[28] );
            // *(gPtr + 58720256) = rPtr_3[7];
            *(gPtr + 58720256) = rPtr[28];
        
            // rPtr_3[8] = *(gPtr + 67108864);
            // turboFFT_ZADD(rPtr_3[8], rPtr_3[8], rPtr[2] );
            // *(gPtr + 67108864) = rPtr_3[8];
            *(gPtr + 67108864) = rPtr[2];
        
            // rPtr_3[9] = *(gPtr + 75497472);
            // turboFFT_ZADD(rPtr_3[9], rPtr_3[9], rPtr[18] );
            // *(gPtr + 75497472) = rPtr_3[9];
            *(gPtr + 75497472) = rPtr[18];
        
            // rPtr_3[10] = *(gPtr + 83886080);
            // turboFFT_ZADD(rPtr_3[10], rPtr_3[10], rPtr[10] );
            // *(gPtr + 83886080) = rPtr_3[10];
            *(gPtr + 83886080) = rPtr[10];
        
            // rPtr_3[11] = *(gPtr + 92274688);
            // turboFFT_ZADD(rPtr_3[11], rPtr_3[11], rPtr[26] );
            // *(gPtr + 92274688) = rPtr_3[11];
            *(gPtr + 92274688) = rPtr[26];
        
            // rPtr_3[12] = *(gPtr + 100663296);
            // turboFFT_ZADD(rPtr_3[12], rPtr_3[12], rPtr[6] );
            // *(gPtr + 100663296) = rPtr_3[12];
            *(gPtr + 100663296) = rPtr[6];
        
            // rPtr_3[13] = *(gPtr + 109051904);
            // turboFFT_ZADD(rPtr_3[13], rPtr_3[13], rPtr[22] );
            // *(gPtr + 109051904) = rPtr_3[13];
            *(gPtr + 109051904) = rPtr[22];
        
            // rPtr_3[14] = *(gPtr + 117440512);
            // turboFFT_ZADD(rPtr_3[14], rPtr_3[14], rPtr[14] );
            // *(gPtr + 117440512) = rPtr_3[14];
            *(gPtr + 117440512) = rPtr[14];
        
            // rPtr_3[15] = *(gPtr + 125829120);
            // turboFFT_ZADD(rPtr_3[15], rPtr_3[15], rPtr[30] );
            // *(gPtr + 125829120) = rPtr_3[15];
            *(gPtr + 125829120) = rPtr[30];
        
            // rPtr_3[16] = *(gPtr + 134217728);
            // turboFFT_ZADD(rPtr_3[16], rPtr_3[16], rPtr[1] );
            // *(gPtr + 134217728) = rPtr_3[16];
            *(gPtr + 134217728) = rPtr[1];
        
            // rPtr_3[17] = *(gPtr + 142606336);
            // turboFFT_ZADD(rPtr_3[17], rPtr_3[17], rPtr[17] );
            // *(gPtr + 142606336) = rPtr_3[17];
            *(gPtr + 142606336) = rPtr[17];
        
            // rPtr_3[18] = *(gPtr + 150994944);
            // turboFFT_ZADD(rPtr_3[18], rPtr_3[18], rPtr[9] );
            // *(gPtr + 150994944) = rPtr_3[18];
            *(gPtr + 150994944) = rPtr[9];
        
            // rPtr_3[19] = *(gPtr + 159383552);
            // turboFFT_ZADD(rPtr_3[19], rPtr_3[19], rPtr[25] );
            // *(gPtr + 159383552) = rPtr_3[19];
            *(gPtr + 159383552) = rPtr[25];
        
            // rPtr_3[20] = *(gPtr + 167772160);
            // turboFFT_ZADD(rPtr_3[20], rPtr_3[20], rPtr[5] );
            // *(gPtr + 167772160) = rPtr_3[20];
            *(gPtr + 167772160) = rPtr[5];
        
            // rPtr_3[21] = *(gPtr + 176160768);
            // turboFFT_ZADD(rPtr_3[21], rPtr_3[21], rPtr[21] );
            // *(gPtr + 176160768) = rPtr_3[21];
            *(gPtr + 176160768) = rPtr[21];
        
            // rPtr_3[22] = *(gPtr + 184549376);
            // turboFFT_ZADD(rPtr_3[22], rPtr_3[22], rPtr[13] );
            // *(gPtr + 184549376) = rPtr_3[22];
            *(gPtr + 184549376) = rPtr[13];
        
            // rPtr_3[23] = *(gPtr + 192937984);
            // turboFFT_ZADD(rPtr_3[23], rPtr_3[23], rPtr[29] );
            // *(gPtr + 192937984) = rPtr_3[23];
            *(gPtr + 192937984) = rPtr[29];
        
            // rPtr_3[24] = *(gPtr + 201326592);
            // turboFFT_ZADD(rPtr_3[24], rPtr_3[24], rPtr[3] );
            // *(gPtr + 201326592) = rPtr_3[24];
            *(gPtr + 201326592) = rPtr[3];
        
            // rPtr_3[25] = *(gPtr + 209715200);
            // turboFFT_ZADD(rPtr_3[25], rPtr_3[25], rPtr[19] );
            // *(gPtr + 209715200) = rPtr_3[25];
            *(gPtr + 209715200) = rPtr[19];
        
            // rPtr_3[26] = *(gPtr + 218103808);
            // turboFFT_ZADD(rPtr_3[26], rPtr_3[26], rPtr[11] );
            // *(gPtr + 218103808) = rPtr_3[26];
            *(gPtr + 218103808) = rPtr[11];
        
            // rPtr_3[27] = *(gPtr + 226492416);
            // turboFFT_ZADD(rPtr_3[27], rPtr_3[27], rPtr[27] );
            // *(gPtr + 226492416) = rPtr_3[27];
            *(gPtr + 226492416) = rPtr[27];
        
            // rPtr_3[28] = *(gPtr + 234881024);
            // turboFFT_ZADD(rPtr_3[28], rPtr_3[28], rPtr[7] );
            // *(gPtr + 234881024) = rPtr_3[28];
            *(gPtr + 234881024) = rPtr[7];
        
            // rPtr_3[29] = *(gPtr + 243269632);
            // turboFFT_ZADD(rPtr_3[29], rPtr_3[29], rPtr[23] );
            // *(gPtr + 243269632) = rPtr_3[29];
            *(gPtr + 243269632) = rPtr[23];
        
            // rPtr_3[30] = *(gPtr + 251658240);
            // turboFFT_ZADD(rPtr_3[30], rPtr_3[30], rPtr[15] );
            // *(gPtr + 251658240) = rPtr_3[30];
            *(gPtr + 251658240) = rPtr[15];
        
            // rPtr_3[31] = *(gPtr + 260046848);
            // turboFFT_ZADD(rPtr_3[31], rPtr_3[31], rPtr[31] );
            // *(gPtr + 260046848) = rPtr_3[31];
            *(gPtr + 260046848) = rPtr[31];
        }}