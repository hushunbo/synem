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

#graph_title = 'Jaccard index of 31 anatomical regions averaged over 612 T1-T2 registrations'
#data_fname = 'data_multi_seg.txt'
#regions_fname = 'xlabels_multi_seg.txt'
#series_fname = 'series_multi_seg.txt'

graph_title = 'Jaccard index of 31 anatomical regions averaged over 306 T1-B0-T1 registrations'
data_fname = 'data.txt'
regions_fname = 'xlabels.txt'
series_fname = 'series.txt'

#data_fname = 'data_multi_seg_NO_ECC.txt'
#regions_fname = 'xlabels_multi_seg_NO_ECC.txt'
#series_fname = 'series_multi_seg_NO_ECC.txt'


#graph_title = ''
#graph_title = 'LPBA40, 56 regions, 1560 registrations'
#data_fname = 'data_mono_short_t_overlap_lpba40.txt'
#regions_fname = 'xlabels_mono_short_t_overlap_lpba40.txt'
#series_fname = 'series_mono_short_t_overlap_lpba40.txt'


data = np.loadtxt(data_fname)
if len(data.shape)==1:
    data = data.reshape((data.shape[0],1))
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
    line, = ax.plot(range(1, nrows+1), data[:,s], linestyle=linestyles[s], marker=markers[s], linewidth=2)
    line.set_label(series[s])
ax.legend(loc=3, fontsize=14)
plt.xticks(np.array(range(1, 1 + nrows)), xlabels, rotation='vertical', horizontalalignment='left', fontsize=18)
plt.xlim(0.75, nrows+0.25)
plt.grid()
plt.ylim(0.05, 0.8)
plt.tight_layout()
plt.ylabel('Jaccard index', fontsize=24)
plt.title(graph_title, fontsize=26)
