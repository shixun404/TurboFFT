
from math import *
import numpy as np
M_PI = 3.141592653589793
def ft_2D_fft_code_gen_upload1(N, N1, N2, num_block, num_thread,
                       radix=2, signal_per_thread=8, if_abft=False):
    exponent = int(log(N, radix))
    N__ = N1
    print(f"########################### N = 2 ** {exponent} ###########################################################")
    print(f"N={N1}, radix={radix}, N1 / radix = {N1/radix}, signal_per_thread={signal_per_thread}")
    plan = []
    twiddle_type = []
    blockdim_x = int(num_thread // (N1 // signal_per_thread))
    blockdim_y = int(N1 // signal_per_thread)
    i = 1
    while i <= N1 / radix:
        if i > N1 / radix:
            break
        elif i >= N1 / signal_per_thread:
            twiddle_type.append(1)
        else:
            twiddle_type.append(0)
        n = 0
        output = f""
        while i < int(N1) and n < int(log(signal_per_thread, radix)):
            output += f"{i}->"
            i *= radix
            n += 1
        print(output)
        plan.append(n)
    print("plan: ", plan)
    print("twid: ",twiddle_type)
    order = []
    for i in range(signal_per_thread * 2):
        order.append(i)
    offset = signal_per_thread
    ft_fft = f'''extern __shared__ float shared[];
__global__ void __launch_bounds__({num_thread}) fft_radix{radix}_logN{exponent}_1''' + '''(float2* inputs, float2* outputs, float2* r_1) {
'''
    ft_fft += '''
    '''
    for i in range(signal_per_thread):
        ft_fft += f'''float2 temp_{i};
    '''
    ft_fft += f'''
    float2* sdata = (float2*)shared;
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int bx = blockIdx.x;
    int tid = tx + ty * blockDim.x;
    int N = {N1};
    int __id[{signal_per_thread}];
    int batch_id = (blockIdx.x * blockDim.x) / {N2};
    
    float2 tmp;
    float2 tmp_angle, tmp_angle_rot;
    int j;
    int k;
    int tmp_id;
    int n = 1, n_global = 1;
    float2 tmp_angle_bk;
    '''
    
    ft_fft += '''
    float2 r[3];
    r[0].x = 1.0f;
    r[0].y = 0.0f;
    r[1].x = -0.5f;
    r[1].y = -0.8660253882408142f;
    r[2].x = -0.5f;
    r[2].y = 0.8660253882408142f;
    // if(tx == 0 && ty == 0 && bx > 3000) printf("bx %d\\n", bx);
    float2 mem_checksum, mem_checksum_t1;
    '''
    
    
    if 2 * N1 // num_thread >= 1 and (2 * N1 // num_thread <= 4):
        ft_fft += f'''
    #if FT==2
    float{2 * N1 // num_thread} tmp_r;
    tmp_r = *(float{2 * N1 // num_thread}*)(((float*)r_1) + tid * {2 * N1 // num_thread});
    *(float{2 * N1 // num_thread}*)(((float*)sdata) + tid * {2 * N1 // num_thread}) = tmp_r;
    // if(bx == 0)printf("%d, hello\\n", tid);
    #endif
    ''' 
    elif 2 * N1 // num_thread > 4:
        ft_fft += f'''
    #if FT==2
    float4 tmp_r;
    '''
        for i in range(2 * N1 // num_thread // 4):
            ft_fft += f'''
    tmp_r = *(float4*)(((float*)r_1) + tid * {2 * N1 // num_thread} + {i} * 4);
    *(float4*)(((float*)sdata) + tid * {2 * N1 // num_thread} + {i} * 4) = tmp_r;
    // if(bx == 0)printf("%d, hello\\n", tid);
    '''
        ft_fft += '''
    #endif
    '''
    else:
        ft_fft += '''
    #if FT==2
    float tmp_r;
    tmp_r = *(((float*)r_1) + tid * {2 * N1 // num_thread});
    *(((float*)sdata) + tid * {2 * N1 // num_thread}) = tmp_r;
    // if(bx == 0)printf("%d, hello\\n", tid);
    #endif
    '''

    
    n = 1
    n_global = 1
    ft_fft += '''
    '''
    for i in range(signal_per_thread):
        ft_fft += f'''temp_{i} = inputs[(ty + {i} * {(N1 // signal_per_thread)}) * {N2} + (tx + bx * blockDim.x) % {N2} + batch_id * {N}];
    '''
    ft_fft += '''
    '''
    for i in range(signal_per_thread):
        ft_fft += f'''__id[{i}] = {i * blockdim_y} + ty;
    '''
    
    
    ft_fft += '''
#if FT==2
mem_checksum.x = 0;
mem_checksum.y = 0;
mem_checksum_t1.x = 0;
mem_checksum_t1.y = 0;
__syncthreads();
'''
    for i in range(signal_per_thread):
        ft_fft += f'''
        // if(bx == 0 && tid == 0)printf("%d, %f %f, hello\\n", __id[{i}], sdata[__id[{i}]].x, sdata[__id[{i}]].y);
        mem_checksum.x += sdata[__id[{i}]].x * temp_{i}.x - sdata[__id[{i}]].y * temp_{i}.y;
        mem_checksum.y += sdata[__id[{i}]].y * temp_{i}.x + sdata[__id[{i}]].x * temp_{i}.y;
    '''
    ft_fft += '''
    // __syncthreads();
    // mem_checksum_t1.x = mem_checksum.x; 
    mem_checksum_t1.y = mem_checksum.y + mem_checksum.x; 
    // mem_checksum_t1.x += __shfl_xor_sync(0xffffffff, mem_checksum_t1.x, 16,32);
    // mem_checksum_t1.x += __shfl_xor_sync(0xffffffff, mem_checksum_t1.x, 8, 32);
    // mem_checksum_t1.x += __shfl_xor_sync(0xffffffff, mem_checksum_t1.x, 4, 32);
    // mem_checksum_t1.x += __shfl_xor_sync(0xffffffff, mem_checksum_t1.x, 2, 32);
    // mem_checksum_t1.x += __shfl_xor_sync(0xffffffff, mem_checksum_t1.x, 1, 32);
    
    mem_checksum_t1.y += __shfl_xor_sync(0xffffffff, mem_checksum_t1.y, 16,32);
    mem_checksum_t1.y += __shfl_xor_sync(0xffffffff, mem_checksum_t1.y, 8, 32);
    mem_checksum_t1.y += __shfl_xor_sync(0xffffffff, mem_checksum_t1.y, 4, 32);
    mem_checksum_t1.y += __shfl_xor_sync(0xffffffff, mem_checksum_t1.y, 2, 32);
    mem_checksum_t1.y += __shfl_xor_sync(0xffffffff, mem_checksum_t1.y, 1, 32);
    #endif
    '''
    
    for stage_id in range(len(plan)):
        if True:
            batch_size = 1 if twiddle_type[stage_id] == 0 else n_global // (N1 // signal_per_thread)
            n = 1
            n_global_ = 1
            for j in range(plan[stage_id]):
                i = 0 + signal_per_thread // radix
                ft_fft += f'''
    j = {int(i / ((signal_per_thread) // radix))};
    k = {i // batch_size } % {n_global_};
    MY_ANGLE2COMPLEX((float)(j * k) * {(-2.0 * M_PI / (radix * n_global_))}f, tmp_angle);
    tmp_angle_bk = tmp_angle;
    '''
                for batch in range(batch_size):
                    ft_fft += '''
                    tmp_angle = tmp_angle_bk;
    '''
                    for k in range(max(1, n // 2)):
                        i = k * batch_size + batch + signal_per_thread // radix
                        ft_fft += f'''
        tmp_angle_rot.x = {cos(- M_PI / float(n)) if k != 0 else 1.}f;
        tmp_angle_rot.y = {sin(- M_PI / float(n)) if k != 0 else 0.}f;
        MY_MUL(tmp_angle, tmp_angle_rot, tmp);
        tmp_angle = tmp;
        tmp_angle_rot.x = tmp_angle.y;
        tmp_angle_rot.y = -tmp_angle.x;
        '''
                        for kk in range(signal_per_thread // radix // batch_size // n):
                            i = (kk * n + k) * batch_size + batch + signal_per_thread // radix
                            ft_fft += f'''
        MY_MUL(temp_{order[signal_per_thread - offset + i]}, tmp_angle, tmp);
        temp_{order[signal_per_thread - offset + i]} = tmp;
        ''' * (2 if if_abft else 1)
                            if n // 2 != 0:
                                i += (n // 2) * batch_size
                                ft_fft += f'''
        MY_MUL(temp_{order[signal_per_thread - offset + i]}, tmp_angle_rot, tmp);
        temp_{order[signal_per_thread - offset + i]} = tmp;
        ''' * (2 if if_abft else 1)
                for batch in range(batch_size):
                    for i in range(signal_per_thread // radix // batch_size):
                        tmp_id_left = ((i // n) * 2 * n + (i % n)) * batch_size + batch
                        tmp_id_right = ((i // n) * 2 * n + (i % n) + n) * batch_size + batch
                        ft_fft += f'''
        tmp = temp_{order[i * batch_size + batch + signal_per_thread - offset]};
        MY_ADD(tmp, temp_{order[i * batch_size + batch + signal_per_thread + int(signal_per_thread / 2) - offset]}, temp_{order[i * batch_size + batch + signal_per_thread - offset]});
        MY_SUB(tmp, temp_{order[i * batch_size + batch + signal_per_thread + int(signal_per_thread / 2) - offset]}, temp_{order[i * batch_size + batch + signal_per_thread + int(signal_per_thread / 2) - offset]});
        ''' * (2 if if_abft else 1)
                        ft_fft += f'''
        tmp_id = __id[{order[i * batch_size + batch + signal_per_thread - offset]}];
        tmp_id = (tmp_id / n_global) * 2 * n_global + (tmp_id % n_global);
        __id[{order[i * batch_size + batch + signal_per_thread - offset]}] = tmp_id;
        __id[{order[i * batch_size + batch + signal_per_thread + int(signal_per_thread / 2) - offset]}] = tmp_id + {n_global};
        '''
                        order[tmp_id_left + offset] = order[i * batch_size + batch + signal_per_thread - offset]
                        order[tmp_id_right + offset] = order[i * batch_size + batch + signal_per_thread + int(signal_per_thread / 2) - offset]
                ft_fft += f'''
        n_global *= 2;
        '''
                offset = 0 if  offset > 0 else signal_per_thread
                n *= radix
                n_global *= radix
                n_global_ *= radix

        if twiddle_type[stage_id] == 0:
            ft_fft += f'''
    __syncthreads();
    '''
            for i in range(signal_per_thread):    
                ft_fft += f'''
    MY_ANGLE2COMPLEX((float)(-M_PI * 2 * ((ty) / {int(N1 // N__)}) * {i}) / (float)({N__}), tmp_angle);
    MY_MUL(temp_{order[signal_per_thread - offset + i]}, tmp_angle, tmp);
    temp_{order[signal_per_thread - offset + i]} = tmp;
    '''
        
            for i in range(signal_per_thread):
                ft_fft += f'''
    sdata[tx + {blockdim_x} * __id[{order[signal_per_thread - offset + i]}]] = temp_{order[signal_per_thread - offset + i]};
    ''' if True else f'''
    sdata[(__id[{order[signal_per_thread - offset + i]}] / 16) * 17 + 
    (__id[{order[signal_per_thread - offset + i]}] % 16)] = temp_{order[signal_per_thread - offset + i]};
    '''
            ft_fft += f'''
    __syncthreads();	
    '''
            for i in range(signal_per_thread):
                ft_fft += f'''
    temp_{i} = sdata[tx + {blockdim_x} * ({i * blockdim_y} + ty)];
    __id[{i}] = ty + {i * blockdim_y};
    ''' if True else f'''
    temp_{i} = sdata[(({i} * blockDim.x + tx) / 16) * 17 +
                        (({i} * blockDim.x + tx) % 16)];
    __id[{i}] = tx + {i} * {N // signal_per_thread};
    '''
                order[i] = i
            offset = signal_per_thread
            N__ = N__ / (2 ** plan[stage_id])
        
        else:    
            offset = 0 if  offset > 0 else signal_per_thread
            ft_fft += '''
            #if FT==2
            mem_checksum.x = 0;
            mem_checksum.y = 0;
            int r_id;
    '''
            for i in range(signal_per_thread):
                temp_id = order[signal_per_thread - offset + i]
                ft_fft += f'''
            r_id = __id[{order[signal_per_thread - offset + i]}] % 3;
            mem_checksum.x += temp_{temp_id}.x * r[r_id].x - temp_{temp_id}.y * r[r_id].y;
            mem_checksum.y += temp_{temp_id}.y * r[r_id].x + temp_{temp_id}.x * r[r_id].y;
    '''
    
            ft_fft += '''
            mem_checksum.y = mem_checksum.y + mem_checksum.x;
            mem_checksum.y += __shfl_xor_sync(0xffffffff, mem_checksum.y, 16, 32);
            mem_checksum.y += __shfl_xor_sync(0xffffffff, mem_checksum.y, 8, 32);
            mem_checksum.y += __shfl_xor_sync(0xffffffff, mem_checksum.y, 4, 32);
            mem_checksum.y += __shfl_xor_sync(0xffffffff, mem_checksum.y, 2, 32);
            mem_checksum.y += __shfl_xor_sync(0xffffffff, mem_checksum.y, 1, 32);
            
            // mem_checksum.x += __shfl_xor_sync(0xffffffff, mem_checksum.x, 16, 32);
            // mem_checksum.x += __shfl_xor_sync(0xffffffff, mem_checksum.x, 8, 32);
            // mem_checksum.x += __shfl_xor_sync(0xffffffff, mem_checksum.x, 4, 32);
            // mem_checksum.x += __shfl_xor_sync(0xffffffff, mem_checksum.x, 2, 32);
            // mem_checksum.x += __shfl_xor_sync(0xffffffff, mem_checksum.x, 1, 32);
            if(tid % 32 == 0){
                mem_checksum.x  = mem_checksum.y;
                mem_checksum.y = mem_checksum.y - mem_checksum_t1.y;
                sdata[tid / 32] = mem_checksum;
            }
            __syncthreads();
            mem_checksum.x = 0;
            mem_checksum.y = 0;
            '''
            
            ft_fft += f'''
            if(tid < {num_thread // 32})
            '''
            ft_fft += f'''
            mem_checksum = sdata[tid];
            '''
            i = num_thread // 32
            while( i > 1):
                i //= 2
                ft_fft += f'''
                // mem_checksum.x += __shfl_xor_sync(0xffffffff, mem_checksum.x, {i}, 32);
                mem_checksum.y += __shfl_xor_sync(0xffffffff, mem_checksum.y, {i}, 32);
        '''
            ft_fft += '''
            // if(mem_checksum.y > 1)printf("%f, %f, %f\\n", temp_0.x, temp_0.y, mem_checksum.y );
            // if(tid == 0 && bx < 128)printf("up1 %f, %f, %f\\n", mem_checksum.x, mem_checksum.y,mem_checksum.y * mem_checksum.y / mem_checksum.x);
        '''
            
            ft_fft += f'''
            temp_0.x += 0.1f * (mem_checksum.x);
            temp_0.y += 0.1f * (mem_checksum.y);
            
            #endif
            #if defined(LOG_ON)
            if(tid == 0 && bx < 128)printf("up1 %f, %f, %f\\n", mem_checksum.x, mem_checksum.y, mem_checksum.y / mem_checksum.x);
            #endif
            '''
                
            
            
            
            for i in range(signal_per_thread):
                ft_fft += f'''
    MY_ANGLE2COMPLEX((float)(-M_PI * 2 * ((tx + bx * {blockdim_x}) % {N2}) * (__id[{order[signal_per_thread - offset + i]}])) / (float)({N}), tmp_angle);
    MY_MUL(temp_{order[signal_per_thread - offset + i]}, tmp_angle, tmp);
    temp_{order[signal_per_thread - offset + i]} = tmp;
    outputs[((tx + bx * {blockdim_x}) % {N2}) + {N2} * __id[{order[signal_per_thread - offset + i]}] + batch_id * {N}] = temp_{order[signal_per_thread - offset + i]};
    '''
            ft_fft += '''
    }
'''
            break
    return ft_fft
