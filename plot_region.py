import matplotlib.pyplot as plt
import matplotlib.cm as cm

import pyregion
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits

# read in the image
cluster_name="a383_mosaic_065mas_acs_wfc3ir_total_drz_20110825.fits"
f_cluster = pyfits.open(cluster_name)

try:
    import pywcsgrid2
    ax=pywcsgrid2.subplot(111, header=f_cluster[0].header)
except ImportError:
    ax=plt.subplot(111)

ax.imshow(f_cluster[0].data, cmap=cm.gray, vmin=0., vmax=0.00038, origin="lower")

reg_name = "abell383_wfc3ir_outline.reg"
r = pyregion.open(reg_name).as_imagecoord(header=f_cluster[0].header)

from pyregion.mpl_helper import properties_func_default

# Use custom function for patch attribute
def fixed_color(shape, saved_attrs):
    
    attr_list, attr_dict = saved_attrs
    attr_dict["color"] = "red"
    kwargs = properties_func_default(shape, (attr_list, attr_dict))
    
    return kwargs

# select region shape with tag=="Group 1"
r1 = pyregion.ShapeList([rr for rr in r if rr.attr[1].get("tag") == "Group 1"])
patch_list1, artist_list1 = r1.get_mpl_patches_texts(fixed_color)

r2 = pyregion.ShapeList([rr for rr in r if rr.attr[1].get("tag") != "Group 1"])
patch_list2, artist_list2 = r2.get_mpl_patches_texts()

for p in patch_list1 + patch_list2:
    ax.add_patch(p)
for t in artist_list1 + artist_list2:
    ax.add_artist(t)
    
plt.savefig("abell383_wfc3ir_outline.png")
plt.show()
