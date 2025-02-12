import seaborn as sns
import matplotlib.pyplot as plt
import torch as th
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
# plt.rcParams["font.family"] = "Times New Roman"
plt.rc('font', size=16,)
cmap = sns.diverging_palette(220, 10, as_cmap=True)
fontsize = 20
logN_beg, logN_end = 1, 26
logBS_beg, logBS_end = 0, 28
cuFFT = th.load('../artifact_data/cuFFT.pt')
TurboFFT = th.load('../artifact_data/TurboFFT.pt')
vkFFT = th.load('../artifact_data/VkFFT_data/vkFFT.pt')

fig, ax = plt.subplots(2, 2, figsize=(8, 8),)
fig.tight_layout()
fig.subplots_adjust(hspace=0, wspace=0.03)
bd = 40
# print(cuFFT.shape)
mask = th.ones(30, 30)
for i in range(26):
        mask[i][:28-i+1] = 0
mask = mask.numpy()
mask = mask[1:26]
label = [
    '(a) TurboFFT',
    '(b) VkFFT',
    '(c) TurboFFT',
    '(d) VkFFT',
]
precision = ['FP32', 'FP64']
cbar_ax = fig.add_axes([.96, .1, .02, .8], label='Overhead vs cuFFT (%)')
for i in range(2):
	heatmap = (TurboFFT[i]/cuFFT[i])[1:26] * 100 - 100
	sns_plot = sns.heatmap(heatmap, cmap=cmap,  mask=mask, vmax=bd, vmin=-bd,center=0,
	# sns_plot = sns.heatmap((TurboFFT[i]/cuFFT[i]) * 100 - 100, cmap=cmap, vmax=50, center=0,
			square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=ax[i, 0],cbar=False)
	# print(vkFFT[:30, :30].shape)
	heatmap = (vkFFT[i]/cuFFT[i])[1:26] * 100 - 100
	sns_plot = sns.heatmap(heatmap, cmap=cmap, mask=mask, vmax=bd, vmin=-bd,center=0,
	# sns_plot = sns.heatmap((vkFFT[i] / cuFFT[i]) * 100 - 100, cmap=cmap, vmax=50, center=0,
			square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=ax[i, 1],
			cbar_ax=cbar_ax if i == 0 else None, cbar=i == 0)
	ax[i, 0].invert_yaxis()
	ax[i, 1].invert_yaxis()
cbar_ax.set_ylabel('Overhead vs cuFFT (%)', fontdict={'fontsize': fontsize, 'weight': 'bold'})
# fig.colorbar(ax.get_children()[0], orientation="horizontal")

for i in range(2):
	b = 0.1
	ax[i, 0].text(.3 + b, .9, label[i * 2], ha='left', va='top', transform=ax[i, 0].transAxes, fontsize=fontsize, fontdict={'weight': 'bold'})
	ax[i, 1].text(.3 + b, .9, label[i * 2 + 1], ha='left', va='top', transform=ax[i, 1].transAxes, fontsize=fontsize, fontdict={'weight': 'bold'})
	ax[i, 0].text(.45 + b, .8, precision[i], ha='left', va='top', transform=ax[i, 0].transAxes, fontsize=fontsize, fontdict={'weight': 'bold'})
	ax[i, 1].text(.425 + b, .8, precision[i], ha='left', va='top', transform=ax[i, 1].transAxes, fontsize=fontsize, fontdict={'weight': 'bold'})
	for j in range(2):
		b = 0.5
		ax[i, j].set_yticks([0 + b, 4 + b, 9+ b, 14+ b, 19+ b, 24+ b])
		ax[i, j].set_yticklabels([1, 5, 10, 15, 20, 25])
		ax[i, j].set_xticks([i + b for i in range(0, 28, 3)])
		ax[i, j].set_xticklabels([i for i in range(0, 28, 3)])
		
		ax[i, j].set_ylabel('log(N)')
		ax[i, j].set_xlabel('log(Batch size)')

plt.savefig('../artifact_figures/figure1.pdf', bbox_inches='tight')
print('./artifact_figures/figure1.pdf')
# plt.savefig('C:/Users/SHIXUNWU/TurboFFT/fig/TurboFFT_vs_VkFFT_vs_cuFFT.png', bbox_inches='tight')