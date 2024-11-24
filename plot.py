import torch as th
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from random import randint
from utils import load_data, load_data_single
import seaborn as sns

# filename = ['../data_0810/benchmark_cufft_A100_fp32.csv', '../data_0810/benchmark_turbofft_A100_fp32.csv']
filenames = [
        # ['../data_0810/benchmark_cufft_A100_fp32.csv', '../data_0810/benchmark_turbofft_A100_fp32.csv'],
        # ['../data_0810/benchmark_cufft_A100_fp64.csv', '../data_0810/benchmark_turbofft_A100_fp64.csv'],
        # ['../data_0810/cuFFT_fp32_cuda_12.0.csv', '../data_0810/turboFFT_fp32_cuda_12.0.csv'],
        ['../data_0810/cuFFT_fp32_cuda_12.0.csv', '../data_0810//turboFFT_fp32_cuda_12.0_40GB_28.csv'],
        ['../data_0810/cuFFT_fp64_cuda_12.0.csv', '../data_0810/turboFFT_fp64_cuda_12.0.csv'],

        # ['../data_0810/cuFFT_fp32_cuda_12.0_40GB_9.csv', '../data_0810/turboFFT_fp32_cuda_12.0_40GB_28.csv'],
]
fig, ax = plt.subplots(2, 2, figsize=(6, 6))
for j in range(2):
    filename = filenames[j]
    cuFFT = th.ones(30, 30)
    TurboFFT = th.ones(30, 30)

    f_cufft = open(filename[0], 'r')
    f_TurboFFT = open(filename[1], 'r')

    lines_cufft = f_cufft.readlines()
    lines_TurboFFT = f_TurboFFT.readlines()

    for i in range(len(lines_cufft)):
        if lines_cufft[i].split(',')[0] != 'cuFFT' and lines_cufft[i].split(',')[0] != 'turboFFT':
        # if lines_TurboFFT[i].split(',')[0] != 'turboFFT':
            continue
        cuFFT[int(lines_cufft[i].split(',')[1])][int(lines_cufft[i].split(',')[2])] = float(lines_cufft[i].split(',')[3])

    # if j == 
    for i in range(len(lines_TurboFFT)):
        if lines_TurboFFT[i].split(',')[0] != 'turboFFT':
            continue
        TurboFFT[int(lines_TurboFFT[i].split(',')[1])][int(lines_TurboFFT[i].split(',')[2])] = float(lines_TurboFFT[i].split(',')[3])
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns_plot = sns.heatmap((TurboFFT/cuFFT) * 100 - 100, cmap=cmap, vmax=20, vmin=-20,center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=ax[j // 2][j % 2],)
    ax[j//2][j % 2].set_title(f'cuda12.{j//2 * 2} FP{32+32*(j%2)}')
