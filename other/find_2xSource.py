from Posterior_analysis.source_pos import *
from Utils.get_res import *
if __name__=="__main__":
    kwres_lns = get_kwres("settings_f160w_CP_logL2_mask31_I_SC.py")["kwargs_results"]["kwargs_lens"]
    sett_2xS = get_setting_module("f160w_CP_logL2_2xSource",1)

    ra_s2,dec_s2 = sett_2xS.source_params[1][0]["center_x"],sett_2xS.source_params[1][0]["center_y"]
    # check that the initial pos is valid -> it is
    #from conversion import  *
    #print(conv_radec_to_xy(sett_2xS,ra_s2,dec_s2))

    sett_true = get_setting_module("f160w_CP_logL2_mask31_I_SC",1)
    print(get_source_gen([ra_s2],[dec_s2],kwres_lns,sett_true))
