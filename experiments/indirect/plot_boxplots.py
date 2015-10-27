cc_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_syncc_scores.txt")
mi_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_synmi_scores.txt")
em_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_synem_scores.txt")
ecc_t1b0 = np.loadtxt("jaccard_boxplot_indirect_t1b0_synecc_scores.txt")
######T1-T1####
scores = [ecc_t1b0, cc_t1b0, em_t1b0, mi_t1b0]
t1t1_names = ["ECC", "CC", "EM", "MI"]
fig, ax1 = plt.subplots(facecolor="white")
boxplot(scores)
ax1.set_axisbelow(True)
#ax1.set_title('Distribution of mean Jaccard scores (averaged over 31 regions per pair). Registration of 306 T1-B0 pairs.', fontsize=20)
ax1.set_ylabel('Mean Jaccard index', fontsize=20)
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_ylim(0.3, 0.8)
xtick_names = setp(ax1, xticklabels=t1t1_names)
setp(xtick_names, rotation=0, fontsize=20)



figure()
idx = 0;
for i in range(3):
    for j in range(3):
        idx += 1
        if scores[i][j] is None:
            continue
        ax = subplot(3, 3, idx)
        boxplot(scores[i][j][-4:])
        ax.set_axisbelow(True)
        ax.set_title(mod_names[i]+" - "+mod_names[j])
        ax.set_ylabel('Mean Jaccard index')
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        ax.set_ylim(0.3, 0.8)
        xtick_names = setp(ax, xticklabels=names[i][j][-4:])
        setp(xtick_names, fontsize=16)


