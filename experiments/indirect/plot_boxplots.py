cc_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_syncc_scores.txt")
mi_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_synmi_scores.txt")
em_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_synem_scores.txt")
ecc_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_synecc_scores.txt")
bfp_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_synbfp_scores.txt")
######T1-T1####
scores = [ecc_t1b0, cc_t1b0, em_t1b0, mi_t1b0, bfp_t1b0]
t1t1_names = ["ECC", "CC", "EM", "MI", "BFP"]
fig, ax1 = plt.subplots(facecolor="white")
ax1.tick_params(axis='both', which='major', labelsize=24)
boxplot(scores)
ax1.set_axisbelow(True)
#ax1.set_title('Distribution of mean Jaccard scores (averaged over 31 regions per pair). Registration of 306 T1-B0 pairs.', fontsize=20)
#ax1.set_ylabel('Mean Jaccard index', fontsize=32)
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_ylim(0.3, 0.8)
xtick_names = setp(ax1, xticklabels=t1t1_names)
setp(xtick_names, rotation=0, fontsize=20)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(48)



