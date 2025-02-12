# TurboFFT tensor(1.1361)
# cuFFT tensor(1.0662)
# TurboFFT tensor(1.1448)
# cuFFT tensor(1.1079)
import seaborn as sns
import matplotlib.pyplot as plt
import torch as th
# plt.rcParams["font.family"] = "Times New Roman"
plt.rc('font', size=16,)
cmap = sns.diverging_palette(220, 10, as_cmap=True)
fontsize = 20
logN_beg, logN_end = 1, 26
logBS_beg, logBS_end = 0, 28
p = 1
precision = ['FP32', 'FP64']
cuFFT = th.load('../artifact_data/cuFFT.pt')[p, 1:26]
TurboFFT = th.load('../artifact_data/TurboFFT.pt')[p, 1:26]
TurboFFT_FT = th.load('../artifact_data/TurboFFT_FT.pt')[p,:, 1:26]
vkFFT = th.load('../artifact_data/VkFFT_data/vkFFT.pt')
nrows = 5
# fig, ax = plt.subplots(nrows, 2, figsize=(8, 2 + nrows * 3),)
# fig.tight_layout()
# fig.subplots_adjust(hspace=0, wspace=0.03)
bd_max = 100
bd_min = -100
print(cuFFT.shape)
mask = th.ones(30, 30)
for i in range(26):
        mask[i][:28-i+1] = 0
mask = mask.numpy()
mask = mask[1:26]
letter_id = 'abcdef'
label = ['TurboFFT Offline', 'TurboFFT w/ FT,\nThread-Level']
for i in range(6):
	label.append(f'({letter_id[i]}) TurboFFT w/ FT,\nthread_BS={2**i}')
TurboFFT_FT_best = th.ones(25, 30) * 1e6
result = []
for i in range(nrows):
	TurboFFT_FT[i][TurboFFT_FT[i] < 0.9* TurboFFT] = 1e6
	tmp_result = []
	heatmap = (th.min(TurboFFT_FT_best, TurboFFT_FT[i])/TurboFFT_FT_best) * 100 - 100
	tmp_result.append(heatmap)
	TurboFFT_FT_best = th.min(TurboFFT_FT_best, TurboFFT_FT[i])
	heatmap = (TurboFFT_FT_best/cuFFT) * 100 - 100
	tmp_result.append(heatmap)
	result.append(tmp_result)
# assert 0
fig = plt.figure(figsize=(15, 10))
fig.tight_layout()
fig.subplots_adjust(hspace=0, wspace=-1)
ncols = 4
gs = fig.add_gridspec(2, ncols, height_ratios=[ncols // 2, 1])
cbar_ax = [None, 
           fig.add_axes([.99, .5, .01, .4], label='Overhead vs cuFFT (%)'),
           None, None, None, 
            fig.add_axes([.99, .1, .01, .17], label='Overhead vs Previous FT (%)'),
            None, None,
            ]
axs = []
axs.append(fig.add_subplot(gs[0, 0:ncols // 2]))  # Merging first three positions for subplot 1
axs.append(fig.add_subplot(gs[0, ncols // 2:ncols]))  # Merging next two positions for subplot 2
label = ['(a) Offline FT-FFT', '(b) TurboFFT w/ FT', '(c) Thread']
letter_id = 'abcdefghi'
for i in range(6):
	label.append(f'({letter_id[i + 3]}) TB-{2**i}')
for i in range(ncols):
    axs.append(fig.add_subplot(gs[1, i]))
id_i = [0, -1, 1, 2, 3, 4, 5, 6]
id_j = [1, 1, 0, 0, 0, 0, 0, 0]
bd_max = [100, 100, 0, 0, 0, 0, 0, 0]
bd = -10
bd_min = [-100, -100, bd, bd, bd, bd, bd, bd]
ftsize = 9
for i in range(ncols // 2 + ncols):
    ax = axs[i]
    heatmap = result[id_i[i]][id_j[i]]
    sns_plot = sns.heatmap(
        heatmap.to(int).clamp(max=300)[:, :28], cmap=cmap, mask=mask[:, :28], vmax=bd_max[i], vmin=bd_min[i],center=0,
        square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=ax,
        cbar_ax=cbar_ax[i], cbar=cbar_ax[i]!=None, annot=i<2, fmt="d", annot_kws={"size": ftsize, "weight": 'bold', 'color': 'black'})
    ax.invert_yaxis()
    title_fontsize = 30 if i < 2 else 30
    if i < 2:
        ax.text(.35, .9, label[i], ha='left', va='top', 
                transform=ax.transAxes, fontsize=title_fontsize, 
                fontdict={'weight': 'bold'})
        b = 0.5
        # ax.set_yticks([i + b for i in range(0, 25)])
        # # ax.set_yticklabels([i + 1 for i in range(0, 25)], fontdict={'fontsize': title_fontsize, 'weight':'bold'}, rotation=0)
        ax.set_yticks([0 + b, 3 + b, 7 + b, 11 + b, 15 + b, 19 + b, 23+ b])
        ax.set_yticklabels([1, 4, 8, 12, 16, 20, 24], fontdict={'fontsize': title_fontsize,})
        ax.set_xticks([0 + b, 4 + b, 8 + b, 12 + b, 16 + b, 20 + b, 24+ b])
        ax.set_xticklabels(['0', '4', '8', '12',  '16', '20', '24'], fontdict={'fontsize': title_fontsize, },)
        ax.text(.45, .8, f'{precision[p]} on A100', ha='left', va='top', transform=ax.transAxes, fontsize=title_fontsize, fontdict={'weight': 'bold'})
    else:
        # ax.set_title(label[i], fontdict={'fontsize': title_fontsize, 'weight': 'bold'})
        ax.text(.3, .92, label[i], ha='left', va='top', 
                transform=ax.transAxes, fontsize=title_fontsize, 
                fontdict={'weight': 'bold'})
        b = 0.5
        ax.set_yticks([0 + b, 7 + b, 15 + b, 23+ b])
        ax.set_yticklabels([1, 8, 16, 24], fontdict={'fontsize': title_fontsize, },)
        ax.set_xticks([0 + b, 8 + b, 16 + b, 24+ b])
        ax.set_xticklabels(['0', '8', '16', '24'], fontdict={'fontsize': title_fontsize, },)
    if i < 3:
        ax.set_ylabel('log(N)', fontsize=title_fontsize)
        
    ax.set_xlabel('log(Batch size)', fontsize=title_fontsize)
    # ax.invert_yaxis()
cbar_ax[1].set_ylabel('Overhead vs cuFFT (%)', fontdict={'fontsize': 30, 'weight': 'bold'})
cbar_ax[ncols // 2 + ncols - 1].set_ylabel('ABFT Stepwise\nOptimization (%)', fontdict={'fontsize': 30, 'weight': 'bold'})
plt.tight_layout()
# fig.savefig(f'../figs/ABFT_stepwise_optimization_A100_FP{32*(p + 1)}.pdf', bbox_inches='tight')
# fig.savefig(f'D:/25_PPOPP_TurboFFT/figs/ABFT_stepwise_optimization_A100_FP{32*(p + 1)}.pdf', bbox_inches='tight')

plt.savefig(f'../artifact_figures/figure17.pdf', bbox_inches='tight')
print('./artifact_figures/figure17.pdf')