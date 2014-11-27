cc_t1t1 = np.loadtxt("jaccard_boxplot_syncc_t1t1.txt")
cc_t1t2 = np.loadtxt("jaccard_boxplot_syncc_t1t2.txt")
cc_t1pd = np.loadtxt("jaccard_boxplot_syncc_t1pd.txt")
cc_t2pd = np.loadtxt("jaccard_boxplot_syncc_t2pd.txt")
ecc_t1t1 = np.loadtxt("jaccard_boxplot_synecc_t1t1.txt")
ecc_t1t2 = np.loadtxt("jaccard_boxplot_synecc_t1t2.txt")
ecc_t1pd = np.loadtxt("jaccard_boxplot_synecc_t1pd.txt")
ecc_t2pd = np.loadtxt("jaccard_boxplot_synecc_t2pd.txt")
em_t1t1 = np.loadtxt("jaccard_boxplot_synem_t1t1.txt")
em_t1t2 = np.loadtxt("jaccard_boxplot_synem_t1t2.txt")
mi_t1t1 = np.loadtxt("jaccard_boxplot_synmi_t1t1.txt")
mi_t1t2 = np.loadtxt("jaccard_boxplot_synmi_t1t2.txt")
mi_t1pd = np.loadtxt("jaccard_boxplot_synmi_t1pd.txt")
mi_t2pd = np.loadtxt("jaccard_boxplot_synmi_t2pd.txt")


######T1-T1####
fig, ax1 = plt.subplots()
boxplot([cc_t1t1, ecc_t1t1, mi_t1t1, em_t1t1])
ax1.set_axisbelow(True)
ax1.set_title('Distribution of mean Jaccard scores (averaged over 31 regions per pair). Registration of 301 T1-T1 pairs.')
ax1.set_ylabel('Mean Jaccard index')
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_ylim(0.3, 0.8)
xtick_names = setp(ax1, xticklabels=["CC", "ECC", "MI", "EM"])
setp(xtick_names, rotation=45, fontsize=16)



######T1-T2####
fig, ax1 = plt.subplots()
boxplot([cc_t1t1, ecc_t1t2, mi_t1t2, em_t1t2, cc_t1t2])
ax1.set_axisbelow(True)
ax1.set_title('Distribution of mean Jaccard scores (averaged over 31 regions per pair). Registration of 612 T1-T2 pairs.')
ax1.set_ylabel('Mean Jaccard index')
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_ylim(0.3, 0.8)
xtick_names = setp(ax1, xticklabels=["Baseline", "ECC", "MI", "EM", "CC"])
setp(xtick_names, rotation=45, fontsize=16)

######T1-PD####
fig, ax1 = plt.subplots()
boxplot([cc_t1t1, ecc_t1pd, mi_t1pd, cc_t1pd])
ax1.set_axisbelow(True)
ax1.set_title('Distribution of mean Jaccard scores (averaged over 31 regions per pair). Registration of 612 T1-PD pairs.')
ax1.set_ylabel('Mean Jaccard index')
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_ylim(0.3, 0.8)
xtick_names = setp(ax1, xticklabels=["Baseline", "ECC", "MI", "CC"])
setp(xtick_names, rotation=45, fontsize=16)

#####T2-PD####
fig, ax1 = plt.subplots()
boxplot([cc_t1t1, ecc_t2pd, mi_t2pd, cc_t2pd])
ax1.set_axisbelow(True)
ax1.set_title('Distribution of mean Jaccard scores (averaged over 31 regions per pair). Registration of 612 T2-PD pairs.')
ax1.set_ylabel('Mean Jaccard index')
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_ylim(0.3, 0.8)
xtick_names = setp(ax1, xticklabels=["Baseline", "ECC", "MI", "CC"])
setp(xtick_names, rotation=45, fontsize=16)

