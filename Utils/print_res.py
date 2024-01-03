import os,dill
import numpy as np
import argparse as ap
from corner import quantile

from Utils.tools import *
from Utils.get_res import * 
from Data.conversion import qphi_from_e1e2
from Utils.statistical_tools import simp, simp_sig

def get_kwres_sigma(setting,cut_mcmc=0,backup_path="backup_results",save=True,read_save=True,convert_ellipt=True):
    filepath = f"{get_savefigpath(setting)}/_kwres_sigma.dll"
    if read_save:
        try:
            kwargs_results_updated = load_whatever(filepath)
            print(f"Read {filepath}.")
            if kwargs_results_updated["convert_ellipt"]!=convert_ellipt:
                raise RuntimeError("Convert_ellipt parameter has changed, rerunning")
            param_mcmc = get_mcmc_prm(setting)
            if convert_ellipt:
                index_e1,index_e2 = param_mcmc.index("e1_lens0"),param_mcmc.index("e2_lens0")
                param_mcmc[index_e1] = "q_lens0"
                param_mcmc[index_e2] = "phi_lens0"
            return kwargs_results_updated,param_mcmc
        except :
            print(f"File {filepath} not found, creating now.")
            pass

    # modified rewrite_read_results.py
    #MCMC sample
    samples_mcmc = get_mcmc_smpl(setting,backup_path)[cut_mcmc:]
    #parameters' name
    param_mcmc   = get_mcmc_prm(setting,backup_path)
    
    kwargs_results_updated={"convert_ellipt":convert_ellipt} 
    n_ra  = 0
    n_dec = 0  
    if convert_ellipt:
        index_e1,index_e2 = param_mcmc.index("e1_lens0"),param_mcmc.index("e2_lens0")
        e1,e2 = samples_mcmc[:,index_e1],samples_mcmc[:,index_e2]
        q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
        samples_mcmc[:,index_e1] = q
        samples_mcmc[:,index_e2] = phi
        param_mcmc[index_e1] = "q_lens0"
        param_mcmc[index_e2] = "phi_lens0"

    for i in range(len(param_mcmc)):
        param_mcmc[index_e2] = "phi_lens0"

    for i in range(len(param_mcmc)):
        if "psi_ext" in param_mcmc[i]:
        # it's by default in radians, have to convert in degrees
            print("WARNING: psi_ext is now in degrees!")
            samples_mcmc[:,i] = samples_mcmc[:,i]*180/np.pi
        min_v,val,max_v = quantile(samples_mcmc[:,i],q=[0.16, 0.5, 0.84])  
        sig_min = val-min_v
        sig_max = max_v-val
        if sig_min<0 or sig_max<0:
            raise RuntimeError(f"sig_min {sig_min} or sig_max {sig_max} are <0, it should be impossible")
        res_i = {"val":val,"sig_min":sig_min,"sig_max":sig_max} 
        if param_mcmc[i]!="ra_image" and param_mcmc[i]!="dec_image":
            kwargs_results_updated[param_mcmc[i]]=res_i
        elif param_mcmc[i]=="ra_image":
            kwargs_results_updated["ra_image_"+str(n_ra)]=res_i
            n_ra+=1
        else:
            kwargs_results_updated["dec_image_"+str(n_dec)]=res_i
            n_dec+=1
    if save:
        with open(filepath,"wb") as f:
            dill.dump(kwargs_results_updated,f)
        print(f"Saved {filepath}")
    return kwargs_results_updated,param_mcmc

def get_Shajib_val(par_i,conv_to_my_FoR=True):
    #'theta_E_lens0', 'gamma_lens0', 'e1_lens0', 'e2_lens0', 'center_x_lens0', 'center_y_lens0',
    #'theta_E_lens1', 'gamma_ext_lens2', 'psi_ext_lens2
    if "theta_E_lens0"==par_i:
        return r"$1.71\pm{0.01}$"
    elif "gamma_lens0"==par_i:
        return r"$1.96\pm0.04$"
    elif "q_lens0"==par_i:
        return r"$0.51\pm0.04$"
    elif "phi_lens0"== par_i:
        # E of N
        val_i = -81
        if conv_to_my_FoR:
            val_i*=-1
            val_i-=90
        return "$"+str(int(np.round(val_i,0)))+r"\pm6$"
    elif  "gamma_ext" in par_i:
        return r"$0.09\pm 0.02$"
    elif "psi_ext" in par_i:
        # E of N
        val_i = -30
        if conv_to_my_FoR:
            val_i*=-1
            val_i-=90
        return "$"+str(int(np.round(val_i,0)))+r"\pm 3$"
    else:
        return ""
        #raise ValueError(f"{par_i} not found")

def get_S22_val(par_i,conv_to_my_FoR=True):
    #'theta_E_lens0', 'gamma_lens0', 'e1_lens0', 'e2_lens0', 'center_x_lens0', 'center_y_lens0',
    #'theta_E_lens1', 'gamma_ext_lens2', 'psi_ext_lens2
    #dra_A,d_dec_A = 1.558,-0.338
    # correct the right image name
    dra_A,d_dec_A = -0.48,1.820
    dra_ml,d_dec_ml = 0.451,-0.218 
    #sig_radec = 6*1e-3 #6 milliarcsec
    
    if "theta_E_lens0"==par_i:
        return r"$1.581_{-0.002}^{+0.003}$"
    elif "gamma_lens0"==par_i:
        return r"$1.92\pm0.03$"#_{-0.03}^{+0.03}$"
    elif "q_lens0"==par_i:
        return r"$0.96\pm0.01$"#_{-0.01}^{+0.01}$"
    elif "phi_lens0"== par_i:
        # N of E
        val_i = -28.1
        #if conv_to_my_FoR:
        #    val_i*=-1
        #    val_i-=90
        return "$"+str(int(np.round(val_i,0)))+r"_{-2.6}^{+4.5}$"
    elif  "gamma_ext" in par_i:
        return r"$0.127\pm0.004$"
    elif "psi_ext" in par_i:
        # N of E
        val_i = -82.4
        #if conv_to_my_FoR:
        #    val_i*=-1
        #    val_i-=90
        return "$"+str(int(np.round(val_i,0)))+r"\pm0.4$"
    elif "center_x_lens0" in par_i:
        #presented as Dra,Ddec wrt to location
        dra = dra_ml-dra_A
        return "$"+str(np.round(dra,3))+r"\pm 0.006$"
    elif "center_y_lens0" in par_i:
        #presented as Dra,Ddec wrt to location
        ddec = d_dec_ml-d_dec_A
        return "$"+str(np.round(ddec,3))+r"\pm 0.006$"
    else:
        return ""
        #raise ValueError(f"{par_i} not found")

def get_unit(prm_i):
    if "q" in prm_i[:2]:
        return r"[\,]"
    elif "phi" in  prm_i or "psi" in prm_i:
        return r"[$^{\circ}$"
    return r"[\'\']"
    
def tau(dt_i,dt_j,sig_i,sig_j): 
    # copied from stnd_handling_data
    #################################
    #simplified version bc considering method A
    tau_val = abs(dt_i-dt_j)/np.sqrt(sig_i**2+sig_j**2) #~ Z val
    return tau_val
    
def get_s22_val(par_i):
    dra_A,d_dec_A = -0.48,1.820
    dra_ml,d_dec_ml = 0.451,-0.218 
    #sig_radec = 6*1e-3 #6 milliarcsec
    
    if "theta_E_lens0"==par_i:
        res_s22 =1.581
        sig_s22 = [0.002,0.003]
    elif "gamma_lens0"==par_i:
        res_s22 =1.92
        sig_s22 = [0.03]
    elif "q_lens0"==par_i:
        res_s22 = 0.96
        sig_s22 = [0.01]
    elif "phi_lens0"== par_i:
        # N of E
        val_i = -28.1
        
        res_s22 = val_i,
        sig_s22 = [2.6,4.5]
    elif  "gamma_ext" in par_i:
        res_s22 = 0.127
        sig_s22 = [0.004]
    elif "psi_ext" in par_i:
        # N of E
        val_i = -82.4
        res_s22 = val_i
        sig_s22 = 0.4
    elif "center_x_lens0" in par_i:
        #presented as Dra,Ddec wrt to location
        dra = dra_ml-dra_A
        res_s22 = dra
        sig_s22 = 0.006
    elif "center_y_lens0" in par_i:
        #presented as Dra,Ddec wrt to location
        ddec = d_dec_ml-d_dec_A
        res_s22 = ddec
        sig_s22 = 0.006
    else:
        return None,None
    return res_s22,sig_s22
    
def get_tension_S22(res_mine,par_i):
      
    res_s22,sig_s22= get_s22_val(par_i)
    if res_s22 is None and sig_s22 is None:
        return ""
    sig_mine = np.mean([res_mine["sig_min"],res_mine["sig_max"]])
    return str(np.round(tau(res_s22,res_mine["val"],np.mean(sig_s22),sig_mine),2))




@check_setting
def print_table_res(setting,cut_mcmc=0,backup_path="backup_results",row_spacing_mm=1.5,
                    convert_ellipt=True,read_save=True,
                    compare_w_Shajib=False,
                    compare_w_S22=False,
                    print_tension_S22=False):
    kwres,param_mcmc = get_kwres_sigma(setting=setting,cut_mcmc=cut_mcmc,backup_path=backup_path,convert_ellipt=convert_ellipt,read_save=read_save)
    if kwres["convert_ellipt"]!=convert_ellipt:
        raise RuntimeError("Inconsistency between convert_ellipt parameter")
    kwres.pop("convert_ellipt")

    newline = "\n"
    print(newline)
    print(newline)
    str_tbl  = r"\begin{table}"
    str_tbl += newline
    #str_tbl += r"\resizebox{\columnwidth}{!}{"
    #str_tbl += newline
    columns  = "llr"
    if compare_w_Shajib or compare_w_S22:
        columns = "llcr"
        if compare_w_S22 and compare_w_Shajib:
            if print_tension_S22:
                columns = "llcccr"
            else:
                columns = "llccr"
    str_tbl +=r"\begin{tabular}{"+columns+r"}"
    str_tbl += newline
    str_tbl += r"\hline"
    str_tbl += newline 
    for i,kwres_i in enumerate(kwres):
        if "light" in param_mcmc[i] or "image" in param_mcmc[i]:
            continue
        res_i   = kwres[kwres_i]
        str_res = simp(res_i["val"],res_i["sig_min"],res_i["sig_max"],ignore_debug=True) # caution!
        if compare_w_Shajib:
            print(param_mcmc[i],get_Shajib_val(param_mcmc[i]))
            str_res+= " & "+get_Shajib_val(param_mcmc[i])
        elif compare_w_S22:
            str_res+= " & "+get_S22_val(param_mcmc[i])
            if print_tension_S22:
                str_res+=" & " +get_tension_S22(res_i,param_mcmc[i])
        if convert_ellipt:
            if param_mcmc[i]=="q_lens0":
                name_prm = r"$q^{\mathrm{ML}}$"
                udm_prm  = "[]" 
            elif param_mcmc[i]=="phi_lens0":
                name_prm = r"$\phi^{\mathrm{ML}}$"
                udm_prm  = r"[$^\circ$]"
            else:
                name_prm,udm_prm = setting.str_param(param_mcmc[i])
        else:
            name_prm,udm_prm = setting.str_param(param_mcmc[i])
        str_tbl += name_prm+" & "+udm_prm+" & "+str_res+r"\\"
        if row_spacing_mm:
            str_tbl+=f"[{row_spacing_mm}mm]"
        str_tbl +=newline
    str_tbl += r"\hline"
    str_tbl +=newline
    str_tbl +=r"\end{tabular} " #}"
    str_tbl +=newline
    str_tbl +=r"\caption{ }"
    str_tbl +=newline
    str_tbl +=r"\end{table}"
    print(str_tbl)
    print(newline)
    print(newline)


if __name__ == '__main__':
        
    present_program(sys.argv[0])
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Simple print of the results given a single setting file")
    parser.add_argument("-cS","--compare_with_Shajib", action="store_true", dest="compare_w_Shajib", default=False,
                    help="Also print in comparison with results from Shajib et al. 2019")
    parser.add_argument("-TS22","--tension_with_s22", action="store_true", dest="tension_S22", default=False,
                    help="Also print the tension with results from STRIDES 2022")
    parser.add_argument("-cS22","--compare_with_S22", action="store_true", dest="compare_w_S22", default=False,
                    help="Also print in comparison with results from STRIDES 2022")
    
    parser.add_argument("-lpr","--load_prev_res", action="store_true", dest="load_prev_res", default=False,
                    help="Load previously computed results (def:false)")
    parser.add_argument("SETTING")
    args = parser.parse_args()
    print_table_res(args.SETTING,read_save=args.load_prev_res,compare_w_Shajib=args.compare_w_Shajib,compare_w_S22=args.compare_w_S22,print_tension_S22=args.tension_S22)
    success(sys.argv[0])
