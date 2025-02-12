import seaborn as sns
import matplotlib.pyplot as plt
import torch as th
import numpy as np
import matplotlib as mpl
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
# mpl.rcParams['text.usetex'] = True

# # Specify the font: Times New Roman
# mpl.rcParams['text.latex.preamble'] = r'\usepackage{times}'

font = {
        'weight' : 'bold',
        'size'   : 20}
mpl.rc('font', **font)
gpu = 'A100'
fontsize = 20
logN_beg, logN_end = 1, 26
logBS_beg, logBS_end = 0, 4
precision = 0
cuFFT = th.load('../artifact_data/cuFFT.pt')[precision]
TurboFFT = th.load('../artifact_data/TurboFFT.pt')[precision]
vkFFT = th.load('../artifact_data/VkFFT_data/vkFFT.pt')[precision]

N_list = [4, 10, 14, 20]
# fig, ax = plt.subplots(2, 3, figsize=(15, 10))
fig = plt.figure(figsize=(15, 10))
ax = []
ax.append(fig.add_subplot(2, 1, 1))
tmp_axs = []
tmp_axs.append(fig.add_subplot(2, 3, 4))
tmp_axs.append(fig.add_subplot(2, 3, 5))
tmp_axs.append(fig.add_subplot(2, 3, 6))
ax.append(tmp_axs)

fig.subplots_adjust(hspace=0.1, wspace=0.15)
color = sns.color_palette(n_colors=5)
# A100_fp32_combine = th.cat([cuFFT[0], TurboFFT[0], vkFFT[0]], dim=0)
width = 0.23
# for tmp_ax in ax[0, 1:]:
#     tmp_ax.axis('off')
def select_data(N, begin, end):
    bars = []
    bar =  cuFFT / TurboFFT * 100
    bars.append(bar)
    turboFFT_bar = th.cat([bar[N[0]][begin:end], bar[N[1]][begin:end], bar[N[2]][begin:end]], dim=0)

    bar =  cuFFT / vkFFT * 100
    bars.append(bar)
    vkFFT_bar = th.cat([bar[N[0]][begin:end], bar[N[1]][begin:end], bar[N[2]][begin:end]], dim=0)
    
    bar = cuFFT / cuFFT * 100
    bars.append(bar)
    cuFFT_bar = th.cat([bar[N[0]][begin:end], bar[N[1]][begin:end], bar[N[2]][begin:end]], dim=0)
    x = np.arange(logN_end - logN_beg)
    n = logN_end - logN_beg
    # ax[0].set_ylabel()
    tmp_ax = ax[0]
    tmp_ax.bar(x - width, width=width, height=bars[1][logN_beg:logN_end, 0], color=color[2], label = 'VkFFT', edgecolor='black')
    tmp_ax.bar(x, width=width, height=bars[0][logN_beg:logN_end, 0], color=color[0], label = 'TurboFFT', edgecolor='black')
    tmp_ax.bar(x + width, width=width, height=bars[2][logN_beg:logN_end, 0], color=color[1], label = 'cuFFT', edgecolor='black')
    tmp_ax.set_xticks(x)
    tmp_ax.set_xticklabels(['$2^{' + f'{i}' +'}$' for i in range(1, n + 1, 1)], fontdict={'fontsize': 16, 'weight':'bold'}, rotation=0)
    tmp_ax.text(.05, .9, f'FP{32*(precision+1)} on {gpu}, higher is better.', ha='left', va='top',
                 transform=tmp_ax.transAxes, fontdict={'weight': 'bold', 'fontsize':22} )
    tmp_ax.text(.5, .2, f'(a) Changing N, fix batch size=1', ha='center', va='top',
                transform=tmp_ax.transAxes, fontdict={'weight': 'bold', 'fontsize':22},
                bbox=dict(facecolor='white', alpha=0.7) )
    tmp_ax.grid(axis='y')
    tmp_ax.set_xlim([-0.5, 24.5])

    x = np.arange((end - begin))
    n = end - begin
    id = 'bcde'
    for i in range(3):
        tmp_ax = ax[1][i]
        tmp_ax.bar(x - width, width=width, height=vkFFT_bar[i * n:(i + 1) * n], color=color[2], label = 'VkFFT', edgecolor='black')
        tmp_ax.bar(x, width=width, height=turboFFT_bar[i * n:(i + 1) * n], color=color[0], label='TurboFFT', edgecolor='black')
        tmp_ax.bar(x + width, width=width, height=cuFFT_bar[i * n:(i + 1) * n], color=color[1], label = 'cuFFT', edgecolor='black')
        tmp_ax.text(.5, .2, f'({id[i]}) Changing batch size,\nfix N=2^' + f'{N[i]}' + '', ha='center', va='top', transform=tmp_ax.transAxes, 
                    fontdict={'weight': 'bold', 'fontsize':22},
                    bbox=dict(facecolor='white', alpha=0.7))
       
        tmp_ax.set_xticks(x)
        tmp_ax.set_xticklabels([f'{2**i}' for i in range(n)], fontdict={'fontsize': 20, 'weight':'bold'}, rotation=30)
        # tmp_ax.set_xlabel(f'Batch size\nFix N = 2^{N[i]}', fontdict={'fontsize': 22, 'weight':'bold'})
        tmp_ax.locator_params(axis='y', nbins=4) 
        tmp_ax.tick_params(axis='y', labelsize=20, labelrotation=90)
        tmp_ax.yaxis.labelpad = 0
        tmp_ax.grid(axis='y')

select_data([10, 20, 14], 0, 8)
ax[1][0].set_ylabel('Performance vs cuFFT (%)', fontdict={'fontsize': 22, 'weight':'bold'})
ax[0].set_ylabel('Performance vs cuFFT (%)', fontdict={'fontsize': 22, 'weight':'bold'})
ax[0].tick_params(axis='y', labelsize=20, labelrotation=90)
if precision == 1:
    ax[0].set_ylim([0, 175])
# ax[1][0].legend( ncol=3, fontsize=22, loc=(0.6, -0.28), frameon=False)
ax[0].legend( ncol=2, fontsize=22,  loc='upper right', frameon=False)
# fig.savefig('../figs/3xbarchart_overhead_A100_FP32.pdf', bbox_inches='tight')
# fig.savefig(f'D:/25_PPOPP_TurboFFT/figs/3xbarchart_overhead_A100_FP{32*(precision + 1)}.pdf', bbox_inches='tight')
plt.savefig(f'../artifact_figures/figure12.pdf', bbox_inches='tight')
print('./artifact_figures/figure12.pdf')
# plt.show()
# for i, N in enumerate(N_list):
#     begin = 0
#     end = 30
#     x = np.arange(begin,end)
#     select_n = N
    
#     print(turboFFT_bar.shape, x.shape)
#     ax[i].bar(x, width=width, height=turboFFT_bar, color=color[0])
#     ax[i].bar(x + width, width=width, height=cuFFT[select_n][begin:end] / cuFFT[select_n][begin:end], color=color[1])
#     ax[i].bar(x - width, width=width, height=cuFFT[select_n][begin:end] / vkFFT[select_n][begin:end], color=color[2])
#     ax[i].set_title(f'N: 2^{select_n}')
