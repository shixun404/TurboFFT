
from math import *
import numpy as np
M_PI = 3.141592653589793
def ft_1D_fft_code_gen(N, blockdim,
                       radix=2, WorkerFFTSize=8, transpose=0, data_type='float', if_abft=False):             
    exponent = int(log(N, radix))
    N__ = N
    print(f"########################### TRANSPOSE N = 2 ** {exponent} ###########################################################")
    print(f"N={N}, radix={radix}, N / radix = {N / radix}, WorkerFFTSize={WorkerFFTSize}")
    plan = []
    NumWorkersPerSignal = int((N // WorkerFFTSize))
    ThreadblockFFTBatchSize = int(blockdim // (N // WorkerFFTSize))
    
    N_tmp = N
    tensor_shape = []
    while N_tmp > 1:
        plan.append(WorkerFFTSize if N_tmp >= WorkerFFTSize else N_tmp)
        tensor_shape.append(plan[-1])
        N_tmp /= WorkerFFTSize
    tensor_shape = tensor_shape.reverse()
        

    order = []
    for i in range(WorkerFFTSize * 2):
        order.append(i)
    offset = WorkerFFTSize

    ft_fft = f'''extern __shared__ {data_type} shared[];
    __global__ void __launch_bounds__({blockdim}) fft_radix{radix}_logN{exponent}''' \
    + f'''({data_type}2* inputs, {data_type}2* outputs, {data_type}2* r_2)''' + ''' {
    '''

    ft_fft += global2reg(tensor_shape, bs, stride_bs, stride_tensor, stride_threadblock, if_shared=False, bs_major=False)

    for i in len(tensor_shape) - 1:
        next_tensor_shape = tensor_reshape(tensor_shape, )
        ft_fft += fft_reg()
        ft_fft += tensor_reshape_shared(tensor_shape, next_tensor_shape)
        ft_fft += shared2reg(next_tensor_shape)
        tensor_shape = next_tensor_shape
    ft_fft += fft_reg()
    ft_fft += reg2global()
    

    fft += epilogue()
    return ft_fft
