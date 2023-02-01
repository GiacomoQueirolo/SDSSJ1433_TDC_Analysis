import copy
import numpy as np

from tools import *
from get_res import *

present_program(sys.argv[0])

sett_nm = "f140w_CP_logL2_FS"
sett = get_setting_module(sett_nm,1)
mcmc = get_mcmc(sett)
mcmc_smpl = mcmc["mcmc_smpl"]
print("mcmc_smpl shape",np.shape(mcmc_smpl))
mcmc_prm  = mcmc["mcmc_prm"]
x_corr, y_corr = sett.x_image[2],sett.y_image[2]

mcmc_corr = copy.copy(mcmc_smpl)
im_arrays =[[],[]]

first_image = None
for i in range(len(mcmc_prm)):
 prm = mcmc_prm[i]
 if "ra" in prm or "dec" in prm or "_x_" in prm or "_y_" in prm or "center" in prm:
    print(prm)
    if "ra" in prm or "_x_" in prm:
        shift = x_corr
    else:
        shift = y_corr
    mcmc_corr.T[i] = mcmc_smpl.T[i]-shift
    if "ra" in prm:
        if first_image is None:
            first_image=i
        im_arrays[0].append(mcmc_corr.T[i])
    elif "dec" in prm:
        im_arrays[1].append(mcmc_corr.T[i])
        
print("mcmc_corr shape",np.shape(mcmc_corr))

im_arrays = np.array(im_arrays) # shape: 2, 4, n / ra-dec,images, steps

print("im_arrays shape",im_arrays.shape)

mean_im = im_arrays.mean(axis=2)
# invert pos of image A and B in the array
ra_ordered  = [0,-0.00221745, -0.76000637, 2.04098104 ] 
dec_ordered = [0,-3.75414086, -2.12910109, -2.17638207 ]

def image_order(ra,dec,ret_order=True):
    if len(ra)!=len(dec):
        raise ValueError("Give array of the same lenght!")
    image_order = []
    for r,d in zip(ra,dec):
        dist =[]
        for ro,do in zip(ra_ordered,dec_ordered):
            dist.append(np.sqrt((r-ro)**2 + (d-do)**2))
        #print(r,d,dist,np.argsort(dist))
        #min_dist = np.argsort(dist)[0]
        min_dist = np.where(dist==min(dist))[0][0]
        image_order.append(min_dist)
    if not np.all(np.unique(image_order,return_counts=True)[1]==np.ones_like(image_order)):
        raise ValueError("Something went wrong, there are multiply images for the same name")
    
    image_order= np.argsort(image_order)
    
    if ret_order: # I only return the order to be applied to the array 
        return image_order
    else:
        return ra[image_order],dec[image_order]


new_order = image_order(*mean_im)
print("new image order",new_order)


im_arrays_corr = [im_arrays[radec][i] for i in new_order for radec in range(2)]
print("im_arrays_corr shape",np.shape(im_arrays_corr))
mcmc_corr.T[first_image:] = im_arrays_corr
new_im  = [{"ra_image":np.array([mcmc_corr.T[mcmc_prm.index('ra_image'):mcmc_prm.index('ra_image')+4]]),
    'dec_image':np.array([mcmc_corr.T[mcmc_prm.index('ra_image')+4:mcmc_prm.index('ra_image')+8]])} ]

print(np.shape(new_im))
from fermat_pot_analysis import gen_mcmc_fermat

mcmc_ferm_corr = gen_mcmc_fermat(mcmc_corr,mcmc_prm,sett)

mcmc_file_name=get_savemcmcpath(sett)+sett_nm.replace(".py","").replace("settings","mcmc_corr_fermat")+".json"
with open(mcmc_file_name, 'w+') as mcmc_file:
    json.dump(mcmc_fermat, mcmc_file)
success(sys.argv[0])