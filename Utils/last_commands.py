import argparse
from datetime import datetime
import os
from Utils.tools import strip_setting_name

# add some final commands
progs=["fermat_pot_analysis",
       "rewrite_read_results",
       "MCMC_posterior",
       "mcmc_behaviour_plot",
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
"""
for dirpath, dirnames, filenames in os.walk("."):
     for fn in filenames:
         if fn.replace(".py","") in progs:
             if dirnames!=[]:
                 print("'"+dirpath+"/"+dirnames[0]+"/"+fn+"'")
             else:
                 print("'"+dirpath+"/"+fn+"'")
"""
progs=['./Plots/Corner_Plots/psos_plot.py'
'./Plots/Corner_Plots/mcmc_behaviour_plot.py'
'./Plots/Corner_Plots/plot_critical_curves.py'
'./Plots/Corner_Plots/t_corner.py'
'./Utils/create_fits.py'
'./Utils/rewrite_read_results.py'
'./MassToLight/pos_img.py'
'./Data/PSF/check_psf.py'
'./Prior/prior_fermat.py'
'./Prior/prior_source.py'
'./Posterior_analysis/PSO/initial_pos_vs_result.py'
'./Posterior_analysis/PSO/MCMC_posterior.py'
'./Posterior_analysis/PSO/mag_remastered.py'
'./Posterior_analysis/PSO/ellipticity_study.py'
'./Posterior_analysis/PSO/fermat_pot_analysis.py'
'./Posterior_analysis/PSO/pso_behaviour.py']



def last_command(setting_name,prog,log=False,run=False):
    setting_name_wo_py = setting_name.replace(".py","")
    str_com = "python "+str(prog).replace(".py","")+".py "+setting_name_wo_py+".py"
    if log:
        now = datetime.now()
        str_com+=" &> logs/"+setting_name_wo_py.replace("settings_","")+"_"+prog+"_"+str(now.strftime("%d%b"))+".log "
    if "fermat_pot"  not in prog and "rewrite_read_results" not in prog :
        str_com+=" &"
    if run:
        os.system(str_com)
    return(str_com)

#bw_ = bash write
"""
def bw_Python_PID(PyProc=""):
    bw_PyPID=f"PyPID=`ps ax | grep 'python {PyProc}' | grep -v grep | awk '{print $1}`"
    return bw_PyPID
"""

def _bw_wait_PyPID(python_prog=""):
    pypid = "`ps ax | grep 'python "+python_prog+"' | grep -v grep | awk '{print $1}'`"
    bw_wait="\nfor pypid in "+pypid+" ; do \n wait $pypid \ndone\n"
    return bw_wait

def bw_wait_PyPID(python_prog):
    if type(python_prog) is list:
        bww = ""
        for pp in python_prog:
            bww+=_bw_wait_PyPID(pp)+"\n"
    elif type(python_prog) is str:
        bww = _bw_wait_PyPID(python_prog)+"\n"
    return bww

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
        
        for ip,prog in enumerate(progs):
            txt+=last_command(setting,prog,log=logs)+"\n"
            if ip==1:
                txt+=bw_wait_PyPID(["fermat_pot_analysis","rewrite_read_results"])
        #log_check=""
        #if logs:
        #    log_check=" logs/"
        txt+="./check_success.py "+stripped+"\n"
        name_txt = "lc_"+stripped+".sh"
        with open(name_txt,"w") as f:
            f.writelines(txt)
        os.system("chmod +x "+name_txt)
        os.system("./"+name_txt+" &")
