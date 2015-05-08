from experiments.registration.dataset_info import *
from experiments.registration.semi_synthetic import *
import nibabel as nib

t1_name = t1_name = get_brainweb("t1", "strip")
t1_nib = nib.load(t1_name)
t1 = t1_nib.get_data().squeeze()

t2_name = t2_name = get_brainweb("t2", "strip")
t2_nib = nib.load(t2_name)
t2 = t2_nib.get_data().squeeze()

m1, v1 = get_mean_transfer(t1, t2)
m2, v2 = get_mean_transfer(t2, t1)

figure()
draw_boxplots(m1, v1)
graph_title = "Estimated T1-T2 transfer. Brainweb template."
plt.title(graph_title, fontsize=40)

figure()
draw_boxplots(m2, v2)
graph_title = "Estimated T2-T1 transfer. Brainweb template."
plt.title(graph_title, fontsize=40)
