#short test for p light/mass
import sys
import argparse
import numpy as np

from tools import *
from get_res import *
from conversion import qphi_from_e1e2
from image_manipulation import extract_q_ll

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Check p component of light and compare it with p mass profile")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings = get_setting_module(args.SETTING_FILES,sett=True)
        
    for sets in settings:
        print(strip_setting_name(sets))
        kw_res = get_kwres(sets)["kwargs_results"]
        e1,e2  = kw_res["kwargs_lens"][0]["e1"],kw_res["kwargs_lens"][0]["e2"]
        q,phi  = qphi_from_e1e2(e1,e2,ret_deg=True)
        if check_if_SUB(sets):
            phi_ll = sets.phi_ll
            try:
                q_ll   = sets.q_ll
            except AttributeError:
                theta_E = sets.lens_params[0][0]["theta_E"]
                q_ll = extract_q_ll(model_name=sets.lens_light_model_name_full,setting=sets,min_ra=theta_E)
        else:
            e1_ll,e2_ll = kw_res["kwargs_lens_light"][0]["e1"],kw_res["kwargs_lens_light"][0]["e2"]
            q_ll,phi_ll = qphi_from_e1e2(e1_ll,e2_ll,ret_deg=True)
        print("q_lens       : ",q)
        print("q_lens_light : ",q_ll)
        diff = q-(q_ll-0.1)
        sigma_q = 0.1
        print("q_lens-(q_lens_light-0.1)",diff)
        if diff>=0:
            print("")
        else:
            print("punishing term : ", -diff**2/sigma_q**2/2)
        
    success(sys.argv[0])