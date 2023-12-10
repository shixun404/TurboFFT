import pandas as pd
import numpy as np
from code_gen_1D_batch import ft_1D_fft_code_gen
if __name__ =="__main__":
    radix = 2
    
    df = pd.read_csv(f'parameter_radix{radix}.csv')
    if_abft = False

    N = 2 **     int(df['logN'][0])
    include_list = '''
        #include "./include/fft.cuh"
    '''
    signal_per_thread = int(df['signal_per_thread_1'][0])
    bx = int(df['blockdim_x_1'][0])
    by = int(df['blockdim_y_1'][0])
    data_type = 'double'
    num_thread = bx * by
    fft_kernel = ft_1D_fft_code_gen(N=N, num_thread=num_thread,
                                radix=2, signal_per_thread=signal_per_thread, transpose=True, data_type=data_type, if_abft=False)
    function_name = f'''ft_fft_radix{radix}_logN{int(df['logN'][0])}_reg{signal_per_thread}_{data_type}'''
    with open(f"../generated/{function_name}.cuh", 'w') as f:
        f.write(fft_kernel)