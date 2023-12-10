
from math import *
import numpy as np
M_PI = 3.141592653589793
def ft_1D_fft_code_gen(N, num_thread,
                       radix=2, signal_per_thread=8, transpose=0, data_type='float', if_abft=False):             
    if True: 
        exponent = int(log(N, radix))
        N__ = N
        print(f"########################### TRANSPOSE N = 2 ** {exponent} ###########################################################")
        print(f"N={N}, radix={radix}, N / radix = {N / radix}, signal_per_thread={signal_per_thread}")
        log_threadx = 0# thread to log
        log_thready = 0# thread to log
        log_block = 0# thread to log
        plan = []
        twiddle_type = []
        blockdim_y = int((N // signal_per_thread))
        blockdim_x = int(num_thread // (N // signal_per_thread))
        i = 1
        while i <= N / radix:
            if i > N / radix:
                break
            elif i >= N / signal_per_thread:
                twiddle_type.append(1)
            else:
                twiddle_type.append(0)
            n = 0
            output = f""
            while i < int(N) and n < int(log(signal_per_thread, radix)):
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
        
#         ft_fft = f'''
# __constant__ float r_1[{2 * N}];
# '''
        ft_fft = f'''extern __shared__ {data_type} shared[];
        
    __global__ void __launch_bounds__({num_thread}) fft_radix{radix}_logN{exponent}''' \
    + f'''({data_type}2* inputs, {data_type}2* outputs, {data_type}2* r_2)''' + ''' {
    '''
        
        ft_fft += f'''
    int tid = threadIdx.x + threadIdx.y * blockDim.x;
    '''
        for i in range(signal_per_thread):
            ft_fft += f'''{data_type}2 temp_{i};
        '''
        ft_fft += f'''
        {data_type}2* sdata = ({data_type}2*)shared;
        int tx = threadIdx.x;
        int ty = threadIdx.y;
        int bx = blockIdx.x;
        int N = {N};
        int __id[{signal_per_thread}];
        {data_type}2 tmp;
        {data_type}2 tmp_angle, tmp_angle_rot;
        int j;
        int k;
        int tmp_id;
        int n = 1, n_global = 1;
        {data_type}2 tmp_angle_bk;
        '''
        n = 1
        n_global = 1
        for i in range(signal_per_thread):
            ft_fft += f'''temp_{i} = inputs[(ty + {i} * {blockdim_y}) + (tx + bx * {blockdim_x}) * {N}];
        '''
        ft_fft += '''
        '''
        for i in range(signal_per_thread):
            ft_fft += f'''__id[{i}] = {i * blockdim_y} + ty;
        '''
        
        
        for stage_id in range(len(plan)):
            if True:
                batch_size = 1 if twiddle_type[stage_id] == 0 else n_global // (N // signal_per_thread)
                n = 1
                n_global_ = 1
                for j in range(plan[stage_id]):
                    i = 0 + signal_per_thread // radix
                    ft_fft += f'''
        j = {int(i / ((signal_per_thread) // radix))};
        k = {i // batch_size } % {n_global_};
        MY_ANGLE2COMPLEX(({data_type})(j * k) * {(-2.0 * M_PI / (radix * n_global_))}f, tmp_angle);
        tmp_angle_bk = tmp_angle;
        ''' * (2 if if_abft else 1)
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
            ''' * (2 if if_abft else 1)
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
                        for i in range(signal_per_thread // radix  // batch_size):
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
    MY_ANGLE2COMPLEX(({data_type})(-M_PI * 2 * ((ty) / {int(N // N__)}) * {i}) / ({data_type})({N__}), tmp_angle);
    MY_MUL(temp_{order[signal_per_thread - offset + i]}, tmp_angle, tmp);
    temp_{order[signal_per_thread - offset + i]} = tmp;
    '''
                for i in range(signal_per_thread):
                    ft_fft += f'''
        sdata[tx + {blockdim_x} * __id[{order[signal_per_thread - offset + i]}]] = temp_{order[signal_per_thread - offset + i]};
        ''' if False else f'''
        sdata[((tx + {blockdim_x} * __id[{order[signal_per_thread - offset + i]}]) / 16) * 17 + 
        ((tx + {blockdim_x} * __id[{order[signal_per_thread - offset + i]}]) % 16)] = temp_{order[signal_per_thread - offset + i]};
        '''
                ft_fft += f'''
        __syncthreads();		
        '''
                for i in range(signal_per_thread):
                    ft_fft += f'''
        temp_{i} = sdata[tx + {blockdim_x} * ({i * blockdim_y} + ty)];
        __id[{i}] = ty + {i * blockdim_y};
        ''' if False else f'''
        temp_{i} = sdata[((tx + {blockdim_x} * ({i * blockdim_y} + ty)) / 16) * 17 +
                            ((tx + {blockdim_x} * ({i * blockdim_y} + ty)) % 16)];
        __id[{i}] = ty + {i * blockdim_y};
        '''
                    order[i] = i
                offset = signal_per_thread
                N__ = N__ / (2 ** plan[stage_id])
            else:
                ft_fft += f'''
        n_global *= 2;
        '''
                for i in range(signal_per_thread):
                    ft_fft += f'''
        outputs[(tx + bx * {blockdim_x}) * {N} +  __id[{order[signal_per_thread - offset + i]}]] = temp_{order[signal_per_thread - offset + i]};
        '''
                ft_fft += '''
        }
    '''
                break
                                
    return ft_fft
