# compare lens light profiles (p,phi and or e1,e2)

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

from tools import *
from get_res import *
from create_fits import *
from conversion import qphi_from_e1e2,e12_from_qphi


from input_data import init_multi_band_list
#### readapted from plotting_tools:
#from plotting_tools import my_corner_general
import scipy.stats as st
from lenstronomy.Plots.model_plot import ModelPlot





if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    #parser = ArgumentParser(description="Plot comparison for ellipticity parameter of lens light results")
    parser = ArgumentParser(description="Plot comparison for ellipticity parameter of lens light results by plotting them")
    parser.add_argument("-n","--name_res",dest="name_res", help="Name of the result png",default="ellipt_comp.png")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting file(s) to model")
    args          = parser.parse_args()
    setting_names = [s.replace(".py","") for s in get_setting_name(args.SETTING_FILES)]
    backup_path   = "backup_results"
    name_res      = backup_path+"/"+str(args.name_res).replace(".png","")+".png"
    """
    kw_qphi = {} #res  = []
    for sett in setting_names:
        setting     = get_setting_module(sett,1)
        print_setting(setting)
        if check_if_SUB(setting):
            #res.append({"dist":False,"setting":strip_setting_name(setting),"q":setting.q_ll,"phi":setting.phi_ll })
            kw_qphi[setting.filter_name]={"q":setting.q_ll,"phi":setting.phi_ll }
        else:
            mcmc_e1,mcmc_e2 = get_mcmc_smpl_for_prm(setting,"e1_lens_light0"),get_mcmc_smpl_for_prm(setting,"e2_lens_light0")
            mcmc_q,mcmc_phi = qphi_from_e1e2(mcmc_e1,mcmc_e2)
            kw_qphi[setting.filter_name]={"q":np.mean(mcmc_q),"phi":np.mean(mcmc_phi) }
            #res.append({"dist":True,"setting":strip_setting_name(setting),"q":mcmc_q,"phi":mcmc_phi})
    samples = [[ri["q"],ri["phi"]] for ri in res if ri["dist"] ]
    s_names = [[ri["setting"]] for ri in res if ri["dist"] ]
    vaxes   = [[{"value":ri["q"],"label":ri["setting"]} for ri in res if not ri["dist"] ],
                [{"value":ri["phi"],"label":ri["setting"]} for ri in res if not ri["dist"] ] ]
    p = remade_corner_general(samples=samples,samples_names=s_names,vaxes=vaxes)
    print("Saving "+name_res)
    p.savefig(name_res)
    """
    f, axes = plt.subplots(3,len(setting_names), figsize=(16, 8), sharex=False, sharey=False)
    kwargs_result = {"kwargs_lens_light":[{"center_x":1,"center_y":-2}]}
    for isett,sett in enumerate(setting_names):
        setting = get_setting_module(sett,1)
        print_setting(setting)
        if check_if_SUB(setting):
            q = setting.q_ll 
            phi = setting.phi_ll
            e1,e2 = e12_from_qphi(phi=phi,q=q)
            R_sersic,n_sersic = .75,6. # placeholders!!
            print("R and n for sersic profiles are set fix with placeholders for the optical images")
        else:
            kw_res_ll = get_kwres(setting)["kwargs_results"]["kwargs_lens_light"][0]
            R_sersic,n_sersic,e1,e2 = kw_res_ll["R_sersic"],kw_res_ll["n_sersic"],kw_res_ll["e1"],kw_res_ll["e2"]
        kwargs_model = {'z_lens':setting.z_lens,
                'z_source':setting.z_source,
                'lens_model_list': None,
                'lens_light_model_list':  ["SERSIC_ELLIPSE"],
                'point_source_model_list': None,
                'additional_images_list': [False], #list of bools (same length as point_source_type_list).
                # If True, search for additional images of the same source is conducted.
                 'fixed_magnification_list': setting.fixed_mag,  # list of bools (same length as point_source_type_list).
                #If True, magnification ratio of point sources is fixed to the one given by the lens model 
                }
        kwargs_result["kwargs_lens_light"][0]["R_sersic"] = R_sersic
        kwargs_result["kwargs_lens_light"][0]["n_sersic"] = n_sersic
        kwargs_result["kwargs_lens_light"][0]["e1"] = e1
        kwargs_result["kwargs_lens_light"][0]["e2"] = e2
        print(kwargs_result)
        multi_band_list,mask = init_multi_band_list(setting=setting,return_mask=True)
        modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,likelihood_mask_list=[mask.tolist()],\
                        arrow_size=0.02, cmap_string="gist_heat")
        modelPlot.decomposition_plot(ax=axes[0,isett], text='Lens light convolved: '+strip_setting_name(setting), lens_light_add=True,
                         v_min=setting.v_min, v_max=setting.v_max)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    print("Saving "+name_res)
    plt.savefig(name_res)

    success(sys.argv[0])
