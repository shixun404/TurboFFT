import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import torch as th
import seaborn as sns
import sys
import os
# sys.path.append(os.path.abspath('../scripts'))
import utils
color = sns.color_palette(n_colors=5)
ms = 9
N = th.as_tensor([i for i in range(25)])
N += 1

plt.rc('font', size=40, weight='bold')
plt.rcParams['lines.linewidth'] = 3
# plt.rcParams["font.family"] = "Times New Roman"
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 7), )
plt.tight_layout()
fig.subplots_adjust(hspace=0.1, wspace = 0.1)


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

_, cufft = utils.load_data_single_sc(file_name_lists[0])
_, turbofft_no_ft = utils.load_data_single_sc(file_name_lists[1])
_, turbofft_ft = utils.load_data_single_sc(file_name_lists[2])
_, turbofft_err = utils.load_data_single_sc(file_name_lists[3])
_, xin = utils.load_data_single_sc(file_name_lists[-2])
# vkfft = th.as_tensor([54.037978,107.172698,162.357973,214.765546,269.936302,325.183609,378.504591,428.645155,481.552961,529.716212,465.167462,318.590441,346.845709,370.003057,395.774621,389.102683,400.783646,360.550859,361.347375,405.507960,355.067390,362.761015,373.276336,362.836422,334.851178])
# vkfft *= 4

vkFFT = th.load('../artifact_data/VkFFT_data/vkFFT.pt')[1]
vkfft = []
for i in range(1, 26):
    gflops = 5 * (2 ** i) * i * (2 ** (28-i)) / vkFFT[i][28 - i] * 1000 / 1000000000.0
    vkfft.append(gflops)
vkfft = th.as_tensor(vkfft)


l = N.shape[0]
# print(((turbofft_err[1:, 0, 1] - turbofft_no_ft[1:, 0, 1]) / turbofft_no_ft[1:, 0, 1]).mean())
# print((turbofft_err[1:, 0, 1] - turbofft_no_ft[1:, 0, 1]))
# assert 0
roofline_model = th.ones_like(N) * 5600
# print_average((abft_baseline - abft_kernel_huge_err_injec - 100) / abft_baseline, "cublas")
# # assert 0
# print_average((kernel_sgemm_huge - abft_kernel_huge_err_injec) / kernel_sgemm_huge, "sgemm")
# print_average((abft_kernel_huge - abft_kernel_huge_err_injec) / abft_kernel_huge, "no error injection")
# print(N.shape, cufft[:, 0, 1].shape)
ax[1].plot(N, cufft[1:, 0, 1] / unit,  label="cuFFT", marker='o', markersize=ms, color = color[1], clip_on=False)
ax[1].plot(N,  vkfft / unit, label="VkFFT", marker='*',markersize=ms,  color =  color[0], clip_on=False)
ax[1].plot(N, turbofft_no_ft[1:, 0, 1] / unit,  label="TurboFFT w/o FT", marker='P',markersize=ms, color = color[2], clip_on=False)
ax[1].plot(N, turbofft_ft[1:, 0, 1] / unit,  label="TurboFFT w/ FT", marker='s',markersize=ms,  color = color[3], clip_on=False)
ax[1].plot(N,  turbofft_err[1:, 0, 1] / unit, '--', label="TurboFFT: err. inj.", marker='^',markersize=ms,  color =  'purple', clip_on=False)
ax[1].plot(N,  xin[1:, 0, 1] / unit, '--', label="Offline FT-FFT: err. inj.", marker='D',markersize=ms,  color =  'k', clip_on=False)


# ax[0].plot(N, roofline_model[:l] / unit,  label="roofline", marker='x', markersize=ms, color = 'k', clip_on=False)


ax[1].set_xlabel("(b) FP64 on A100",fontdict=dict(weight='bold',  size=40))
ax[1].text(.90, -.04, 'logN', ha='left', va='top', transform=ax[1].transAxes)

# ax[0].set_xscale("log", base=2)
# ax[0].set_xlim(0, 10240)
# xticks = [0, 2000, 4000, 6000, 8000, 10000]#[256, 512, 1024, 2048, 4096, 8192]
# ax[0].set_xticklabels(['0','2k','4k','6k','8k','10k',])
# ax[0].set_xticks(xticks)
# yticks = [1, 2, 3, 4, 5, 5.5]
# ax[0].set_yticks(yticks)
# ylabels = ['1', '2', '3', '4', '5', '5.5']
# ax[0].set_yticklabels(ylabels)
# ax[0].set_xticklabels([f'{s}' for s in xticks])

ax[1].grid()

ax[1].legend(loc="lower right", prop={'size': 24, }, labelspacing=0,bbox_to_anchor=(1.02,-0.03))

# roofline_model = th.ones_like(N) * 5600

_, cufft = utils.load_data_single_sc(file_name_lists[0 + 4])
_, turbofft_no_ft = utils.load_data_single_sc(file_name_lists[1 + 4])
_, turbofft_ft = utils.load_data_single_sc(file_name_lists[2 + 4])
_, turbofft_err = utils.load_data_single_sc(file_name_lists[3 + 4])
_, xin = utils.load_data_single_sc(file_name_lists[-1])
# vkfft  = th.as_tensor([106.505101,211.794537,320.018956,426.183860,534.689379,645.786371,744.521505,856.230688,965.247954,1080.585532,1167.387529,1053.020716,680.060203,730.905417,776.663035,825.653568,892.905666,822.305245,946.622101,894.953900,714.654186,748.555643,776.556882,806.125929,756.330899])
# vkfft *= 4

vkFFT = th.load('../artifact_data/VkFFT_data/vkFFT.pt')[0]
vkfft = []
for i in range(1, 26):
    gflops = 5 * (2 ** i) * i * (2 ** (28-i)) / vkFFT[i][28 - i] * 1000 / 1000000000.0
    vkfft.append(gflops)
vkfft = th.as_tensor(vkfft)

l = N.shape[0]
# print(((turbofft_err[1:, 0, 1] - turbofft_no_ft[1:, 0, 1]) / turbofft_no_ft[1:, 0, 1]).mean())
# print((turbofft_err[1:, 0, 1] - turbofft_no_ft[1:, 0, 1]))
# assert 0
roofline_model = th.ones_like(N) * 5600
# print_average((abft_baseline - abft_kernel_huge_err_injec - 100) / abft_baseline, "cublas")
# # assert 0
# print_average((kernel_sgemm_huge - abft_kernel_huge_err_injec) / kernel_sgemm_huge, "sgemm")
# print_average((abft_kernel_huge - abft_kernel_huge_err_injec) / abft_kernel_huge, "no error injection")

ax[0].plot(N, cufft[1:, 0, 1] / unit,  label="cuFFT", marker='o', markersize=ms, color = color[1], clip_on=False)
ax[0].plot(N,  vkfft / unit, label="VkFFT", marker='*',markersize=ms,  color =  color[0], clip_on=False)
ax[0].plot(N, turbofft_no_ft[1:, 0, 1] / unit,  label="TurboFFT w/o FT", marker='P',markersize=ms, color = color[2], clip_on=False)
ax[0].plot(N, turbofft_ft[1:, 0, 1] / unit,  label="TurboFFT w/ FT", marker='s',markersize=ms,  color = color[3], clip_on=False)
ax[0].plot(N,  turbofft_err[1:, 0, 1] / unit, '--', label="TurboFFT: err. inj.", marker='^',markersize=ms,  color =  'purple', clip_on=False)
ax[0].plot(N,  xin[1:, 0, 1] / unit, '--', label="Offline FT-FFT: err. inj.", marker='D',markersize=ms,  color =  'k', clip_on=False)
# # ax[1].plot(N, roofline_model[:l] / unit,  label="roofline", marker='x', markersize=ms, color = 'k', clip_on=False)

ax[0].set_xlabel("(a) FP32 on A100", fontdict=dict(weight='bold', size=40))
ax[0].set_ylabel("Perf. (TFLOPS)", fontdict=dict(weight='bold'))
ax[0].text(.90, -.04, 'logN', ha='left', va='top', transform=ax[0].transAxes)
# ax[0].text(.90, -.14, 'BS=2^28/N', ha='left', va='top', transform=ax[0].transAxes)
# ax[1].set_ylabel("Performance (TFLOPS)", fontdict=dict(weight='bold'))

# ax[1].set_xlim(0, 10240)
# xticks = [0, 2000, 4000, 6000, 8000, 10000]#[256, 512, 1024, 2048, 4096, 8192]
# ax[1].set_xticklabels(['0','2k','4k','6k','8k','10k',])
# ax[1].set_xticks(xticks)
# # ax[1].set_xticklabels([f'{s}' for s in xticks])
# # yticks = [1, 2, 3, 4, 5, 5.6]
# # ax[1].set_yticks(yticks)
# ylabels = ['1', '2', '3', '4', '5', '8.1']
# # ax[1].set_yticklabels(ylabels)

ax[0].grid()
yticks = [ 0.5, 1, 1.5, 2, 2.5]
ax[1].set_yticks(yticks)
yticks = [0, 1, 2, 3, 4, 5]
ax[0].set_yticks(yticks)
xticks = [ 0, 5, 10, 15, 20 ]
ax[1].set_xticks(xticks)
ax[0].set_xticks(xticks)
ax[0].yaxis.set_tick_params(rotation=90)
ax[1].yaxis.set_tick_params(rotation=90)

ax[0].legend(loc="lower right", prop={'size': 24},  labelspacing=0, bbox_to_anchor=(1.02,-0.03))
fig.savefig(f"../artifact_figures/figure20.pdf", bbox_inches='tight')
print('./artifact_figures/figure20.pdf')