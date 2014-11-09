def compute_ranks(T):
    n = T.shape[0]
    m = T.shape[1]
    R = np.zeros_like(T, dtype=np.int32)
    for i in range(n):
        for j in range(m):
            rank = 1
            for k in range(m):
                if(j ==k):
                    continue
                if T[i, k] > T[i, j]:
                    rank += 1
            R[i, j] = rank
    return R

#graph_title = 'Jaccard index of 31 anatomical regions averaged over 612 multi-modal registrations'
graph_title = ''
#graph_title = 'Jaccard index of 31 anatomical regions averaged over 306 mono-modal registrations'
#data_fname = 'data_mono_seg.txt'
#regions_fname = 'xlabels_mono_seg.txt'
#series_fname = 'series_mono_seg.txt'
data_fname = 'data_multi_seg_NO_ECC.txt'
regions_fname = 'xlabels_multi_seg_NO_ECC.txt'
series_fname = 'series_multi_seg_NO_ECC.txt'


data = np.loadtxt(data_fname)
nrows, ncols = data.shape

with open(regions_fname) as input:
    xlabels = [s for s in input.readlines()]
    
with open(series_fname) as input:
    series = input.readlines()

markers = ['o','D','s','^']
linestyles = ['--', '--', '--', '--']

fig = plt.figure()
ax = fig.add_subplot(111)
for s in range(ncols):
    line, = ax.plot(range(1, nrows+1), data[:,s], linestyle=linestyles[s], marker=markers[s])
    line.set_label(series[s])
ax.legend()
plt.xticks(np.array(range(1,32)), xlabels, rotation='vertical', horizontalalignment='left')
plt.grid()
plt.ylim(0, 1)
plt.tight_layout()
plt.title(graph_title, fontsize=32)
plt.ylabel('Jaccard index', fontsize=24)
