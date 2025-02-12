import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import torch as th
import seaborn as sns
import sys
import os
sys.path.append(os.path.abspath('../scripts'))
import utils
color = sns.color_palette(n_colors=5)
ms = 9
N = th.as_tensor([i for i in range(25)])
N += 1

plt.rc('font', size=40, weight='bold')
plt.rcParams['lines.linewidth'] = 3
plt.rcParams["font.family"] = "Times New Roman"
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 8), )
plt.tight_layout()
fig.subplots_adjust(hspace=0.1, wspace = 0.2)


unit = 1000
def print_average(tensor, name):
    avg = (tensor.sum() - tensor.max() - tensor.min()) / tensor.shape[0]
    print(name, avg)
    print(tensor.sort())
file_name_lists = ['../artifact_data/cuFFT_data/Benchmark=2_cuFFT_FP64.csv', 
                  '../artifact_data/TurboFFT_data/Benchmark=2_TurboFFT_FP64.csv',
                  '../artifact_data/TurboFFT_data/Benchmark=2_TurboFFT_FP64_FT_BS=1.csv',
                  '../artifact_data/TurboFFT_data/Benchmark=2_TurboFFT_FP64_ERR_BS=32.csv',
                  '../artifact_data/cuFFT_data/Benchmark=2_cuFFT_FP32.csv', 
                  '../artifact_data/TurboFFT_data/Benchmark=2_TurboFFT_FP32.csv',
                  '../artifact_data/TurboFFT_data/Benchmark=2_TurboFFT_FP32_FT_BS=1.csv',
                  '../artifact_data/TurboFFT_data/Benchmark=2_TurboFFT_FP32_ERR_BS=32.csv',
                  '../artifact_data/cuFFT_data/Benchmark=2_Offline_FP64_ERR.csv',
                  '../artifact_data/cuFFT_data/Benchmark=2_Offline_FP32_ERR.csv',]





_, cufft = utils.load_data_single_sc(file_name_lists[0 + 4])
_, turbofft_no_ft = utils.load_data_single_sc(file_name_lists[1 + 4])
_, turbofft_ft = utils.load_data_single_sc(file_name_lists[2 + 4])
_, turbofft_err = utils.load_data_single_sc(file_name_lists[3 + 4])
_, xin = utils.load_data_single_sc(file_name_lists[-1])


vkFFT = th.load('../artifact_data/VkFFT_data/vkFFT_T4.pt')[0]
vkfft_gflops = []
vkfft_bandwidth = []
for i in range(1, 26):
    # print(vkFFT[i])
    gflops = 5 * (2 ** i) * i * (2 ** (26-i)) / vkFFT[i][26 - i] * 1000 / 1000000000.0
    vkfft_gflops.append(gflops)
    bandwidth = ((2 ** i) * (2 ** (26-i)) * 8 * 2) / vkFFT[i][26 - i] * 1000.0 / 1000000000.0
    vkfft_bandwidth.append(bandwidth)
vkfft = [th.as_tensor(vkfft_gflops), th.as_tensor(vkfft_bandwidth)]


def data_analysis(d1, cufft, name):
    rel = (d1 - cufft) / cufft
    rel_mean = rel.mean()  # Assuming rel.mean() returns a float
    rel_std = rel.std()
    rel_max = rel.max()    # Assuming rel.max() returns a float
    rel_min = rel.min()    # Assuming rel.min() returns a float
    print(f"{name:<20} {rel_mean:<15.2f} {rel_std:<15.2f} {rel_max:<15.2f} {rel_min:<15.2f}")
print_name = ['vkfft',
              'Turbofft no  FT',
              'TurboFFT w/ FT',
              'TurboFFT err',
              'Xin err' ]
data = [
    vkfft[0],
    turbofft_no_ft[1:, 0, 1] ,
    turbofft_ft[1:, 0, 1] ,
    turbofft_err[1:, 0, 1] ,
    xin[1:, 0, 1] 
]
# for i in range(5):
#     data_analysis(data[i], cufft[1:, 0, 1], print_name[i])


l = N.shape[0]
# print(((turbofft_err[1:, 0, 1] - turbofft_no_ft[1:, 0, 1]) / turbofft_no_ft[1:, 0, 1]).mean())
# print((turbofft_err[1:, 0, 1] - turbofft_no_ft[1:, 0, 1]))


l = N.shape[0]

roofline_model = th.ones_like(N) * 5600

x_label = ["(a) Compute Performance", "(b) Memory Throughput"]
y_label = ["TFLOPS", "GB/s"]
units = [1000,1]
for i in range(2):
    unit = units[i]
    ax[i].plot(N, cufft[1:, 0, 1+i] / unit,  label="cuFFT", marker='o', markersize=ms, color = color[1], clip_on=False)
    ax[i].plot(N,  vkfft[i] / unit, label="VkFFT", marker='*',markersize=ms,  color =  color[0], clip_on=False)
    ax[i].plot(N, turbofft_no_ft[1:, 0, 1+i] / unit,  label="TurboFFT w/o FT", marker='P',markersize=ms, color = color[2], clip_on=False)
    ax[i].plot(N, turbofft_ft[1:, 0, 1+i] / unit,  label="TurboFFT w/ FT", marker='s',markersize=ms,  color = color[3], clip_on=False)
    ax[i].plot(N,  turbofft_err[1:, 0, 1+i] / unit, '--', label="TurboFFT: err. inj.", marker='^',markersize=ms,  color =  'purple', clip_on=False)
    ax[i].plot(N,  xin[1:, 0, 1+i] / unit, '--', label="Offline FT-FFT: err. inj.", marker='D',markersize=ms,  color =  'k', clip_on=False)

    ax[i].set_xlabel(x_label[i], fontdict=dict(weight='bold', size=40))
    ax[i].set_ylabel(y_label[i], fontdict=dict(weight='bold'))
    ax[i].text(.90, -.04, 'logN', ha='left', va='top', transform=ax[0].transAxes)

    ax[i].grid()
    yticks = [ 0.5, 1, 1.5, 2, 2.5]
    xticks = [ 0, 5, 10, 15, 20 ]

    ax[i].set_xticks(xticks)


    yticks = [ 1, 2, 3, 4,5]
    # ax[0].set_yticks(yticks)
    ax[i].yaxis.set_tick_params(rotation=90)
ax[0].legend(loc="lower right", prop={'size': 24},  labelspacing=0, bbox_to_anchor=(1.02,-0.03))
fig.savefig(f"../artifact_figures/figure22.pdf", bbox_inches='tight')
print('./artifact_figures/figure22.pdf')