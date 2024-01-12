import pandas as pd
import numpy as np
from code_gen_1D_single_batch_double import ft_1D_fft_code_gen
if __name__ =="__main__":
    radix = 2
    
    df = pd.read_csv(f'parameter_radix{radix}.csv')
    if_abft = False

    N = 2 ** int(df['logN'][0])
    include_list = '''
        #include "./include/fft.cuh"
    '''
    n1 = int(df['signal_per_thread_1'][0])
    bx = int(2 ** int(df['logN1'][0])) // int(df['n1'][0])
    bs = int(df['bs1'][0])
    
    data_type = 'double'
    
    fft_kernel = ft_1D_fft_code_gen(N=N, blockdim=bx,
                                radix=2, WorkerFFTSize=n1, data_type=data_type, if_abft=False)
    function_name = f'''ft_fft_radix{radix}_logN{N}_reg{n1}_{data_type}'''
    with open(f"../generated/{function_name}.cuh", 'w') as f:
        f.write(fft_kernel)