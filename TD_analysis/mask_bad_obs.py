## see Seeing_Transp.ipynb in pycs my pc
import numpy as np
def mask_bad_nights(lcs,delete=True):
    prop  = lcs[0].properties
    mask  = lcs[0].mask
    for i,pr in enumerate(prop):
        if pr["Comment"] != "none":
            mask[i]=False
        
        else:
            if pr["Transp"]!="___":
                if float(pr["Transp"])<.7:
                    mask[i]=False
    for lc in lcs:
        if not delete:
            lc.mask = mask
        else:
            lc.remove_epochs(np.where(mask==False)[0].tolist())
    return lcs