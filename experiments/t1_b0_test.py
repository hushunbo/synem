import nibabel as nib
import experiments.registration.regviz as rv
import experiments.registration.dataset_info as info

base_dir = 'D:/opt/registration/data/SCIL/SCIL_01/warp/'
upt1_fname = 'warpedDiff_SCIL_01_b0_up_strip_SCIL_01_t1_strip.nii.gz'
methods = ['SyNECC', 'SyNCC', 'SyNMI', 'SyNEM']
fnames_upt1 = {}
nibs_upt1 = {}
data_upt1 = {}
for method in methods:
    fnames_upt1[method] = base_dir + method + '/' + upt1_fname
    nibs_upt1[method] = nib.load(fnames_upt1[method])
    data_upt1[method] = nibs_upt1[method].get_data()

t1_name = info.get_scil(1, 't1_strip')
t1_nib = nib.load(t1_name)
t1 = t1_nib.get_data().squeeze()

rv.overlay_slices_with_contours(data_upt1['SyNECC'], t1, slice_type=0, ncontours=7)
rv.overlay_slices_with_contours(data_upt1['SyNECC'], t1, slice_type=1, ncontours=10)
rv.overlay_slices_with_contours(data_upt1['SyNECC'], t1, slice_type=2, slice_index=85, ncontours=6)

rv.overlay_slices_with_contours(data_upt1['SyNMI'], t1, slice_type=0, ncontours=7)
rv.overlay_slices_with_contours(data_upt1['SyNMI'], t1, slice_type=1, ncontours=10)
rv.overlay_slices_with_contours(data_upt1['SyNMI'], t1, slice_type=2, slice_index=85, ncontours=6)

rv.overlay_slices_with_contours(data_upt1['SyNCC'], t1, slice_type=0, ncontours=7)
rv.overlay_slices_with_contours(data_upt1['SyNCC'], t1, slice_type=1, ncontours=10)
rv.overlay_slices_with_contours(data_upt1['SyNCC'], t1, slice_type=2, slice_index=85, ncontours=6)

rv.overlay_slices_with_contours(data_upt1['SyNEM'], t1, slice_type=0, ncontours=7)
rv.overlay_slices_with_contours(data_upt1['SyNEM'], t1, slice_type=1, ncontours=10)
rv.overlay_slices_with_contours(data_upt1['SyNEM'], t1, slice_type=2, slice_index=85, ncontours=6)