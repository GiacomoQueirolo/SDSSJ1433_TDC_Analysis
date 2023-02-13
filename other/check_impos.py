import json
import numpy as np
from corner import corner
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
if __name__=="__main__":

    sett_nm="f140w_CP_logL2"
    with open("mcmc_smpl_"+sett_nm+".json","r") as f:
        mcmc=np.array(json.load(f))

    with open("mcmc_prm_"+sett_nm+".dat","r") as f:
        prm=f.readlines()
    prm= [data_l.replace(",\n","") for data_l in prm]

    index1 = prm.index("ra_image")
    mcmc_pos = mcmc.T[index1:]
    prm_pos  = prm[index1:]
    means    = mcmc_pos.mean(axis=1)
    fg, ax = plt.subplots(len(prm_pos),len(prm_pos),figsize=(15,15))

    other_means = False
    if other_means==True:
        sett_nm2 = "f140w_CP_logL2_FS"
        with open("../mcmc_"+sett_nm2+"/mcmc_smpl_"+sett_nm2+".json","r") as f:
            mcmc2=np.array(json.load(f))

        with open("../mcmc_"+sett_nm2+"/mcmc_prm_"+sett_nm2+".dat","r") as f:
            prm2=f.readlines()
        prm2 = [data_l.replace(",\n","") for data_l in prm2]

        index2 = prm2.index("ra_image")
        mcmc_pos2 = mcmc2.T[index2:]
        means    = mcmc_pos2.mean(axis=1)

    corner(np.array(mcmc_pos).T,labels=prm,show_titles=True,truths=means,truth_color="r",fig=fg)
    means_str = str([np.round(i,2) for i in means]).replace("[","").replace("]","")
    if other_means:
        label="Means of "+sett_nm2
    else:
        label="Means"
    legend_elements = [Patch(label=label,facecolor="r"),Patch(label=means_str,facecolor="w")]
    axdel=ax[0][-1]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    plt.tight_layout()
    plt.savefig("impos.pdf")
    print("Produced impos.pdf")
