import argparse
from datetime import datetime
import os
from tools import strip_setting_name

# add some final commands
progs=["fermat_pot_analysis",
       "MCMC_posterior",
       "mcmc_behaviour_plot",
       "rewrite_read_results",
       "mag_remastered",
       "t_corner",
       "psos_plot",
       "prior_fermat",
       "pso_behaviour",
       "ellipticity_study",
       "pos_img",
       "initial_pos_vs_result",
       "create_fits",
       "plot_critical_curves",
       "prior_source",
       "check_psf"]

def last_comand(setting_name,prog,log=False):
    setting_name_wo_py = setting_name.replace(".py","")
    str_com = "python "+str(prog).replace(".py","")+".py "+setting_name_wo_py+".py"
    if log:
        now = datetime.now()
        str_com+=" &> logs/"+setting_name_wo_py.replace("settings_","")+"_"+prog+"_"+str(now.strftime("%d%b"))+".log "
    if "fermat_pot" not in prog:
        str_com+" &"
    return(str_com)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Run all the complementary programs")
    parser.add_argument("-l","--print_log",action="store_true",dest="logs",default=False,help="Produce log files in the logs/ directory")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings=args.SETTING_FILES
    logs = args.logs
    for setting in settings:
        txt = ""    
        stripped = strip_setting_name(setting)
        for prog in progs:
            txt+=last_comand(setting,prog,log=logs)+"\n"
        
        txt+="./check_success.py "+stripped+"\n"
        name_txt = "lc_"+stripped+".sh"
        with open(name_txt,"w") as f:
            f.writelines(txt)
        os.system("chmod +x "+name_txt)
        os.system("./"+name_txt+" &")