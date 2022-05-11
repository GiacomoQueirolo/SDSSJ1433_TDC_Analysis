import argparse
from datetime import datetime
import os

# add some final commands
progs=["fermat_pot_analysis",
       "MCMC_posterior",
       "mcmc_behaviour_plot",
       "rewrite_read_results",
       "mag_remastered",
       "t_corner",
       "psos_plot",
       "Prior",
       "pso_behaviour",
       "pos_img",
       "ellipticity_study"]

def last_comand(setting_name,prog,log=False):
    setting_name_wo_py = setting_name.replace(".py","")
    str_com = "python "+str(prog).replace(".py","")+".py "+setting_name_wo_py+".py"
    if log:
        now = datetime.now()
        str_com+=" &> logs/"+setting_name_wo_py.replace("settings_","")+"_"+prog+"_"+str(now.strftime("%d%b"))+".log "
    if "fermat_pot" not in prog:
        str_com+" &"
    print(str_com)
    os.system(str_com)
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Run all the complementary programs")
    parser.add_argument("-l","--print_log",action="store_true",dest="logs",default=False,help="Produce log files in the logs/ directory")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings=args.SETTING_FILES
    logs = args.logs
    for setting in settings:    
        for prog in progs:
            last_comand(setting,prog,log=logs)
        os.system("./check_success.py "+setting)
