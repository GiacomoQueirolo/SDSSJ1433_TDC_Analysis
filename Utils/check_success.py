#!/home/moon/queirolo/anaconda3/envs/queirolo_env/bin/python

from argparse import ArgumentParser
import re
import glob
import os
from Utils.tools import *
import sys

def from_log_to_sett(log):
    log_to_sett = log[4:]
    if ".dat" == log_to_sett[-4:]:
        log_to_sett=log_to_sett[:-4]
    log_to_sett = "_".join(log_to_sett.split("_")[:-1])
    log_to_sett = "settings_"+log_to_sett+".py"
    return log_to_sett

def find_latest_log(settings,ignore_log=[]):
    if type(settings)==str:
        settings=[settings]
    files = next(os.walk("."), (None, None, []))[2]
    logs  = [f for f in files if f[:3]=="log" ]
    res_logs =[]
    for sett in settings:
        log = [l for l in logs if strip_setting_name(sett) in l]
        for igl in ignore_log:
            log.remove(igl)    
        if len(log)==1:
            res_logs.append(log[0])
        else:
            mod_time = [os.path.getmtime(l) for l in log]
            sorted_log = [l for _,l in sorted(zip(mod_time,log),key=lambda pair:pair[0])]
            latest_log = sorted_log[-1]
            res_logs.append(latest_log)
    if res_logs==[]:
        raise RuntimeError(f"No logs found for setting(s) {settings}")
    if len(settings)==1:
        return res_logs[0]
    return res_logs

def _check_pso_conv(log,save=True,verbose=False):
    with open(log,"r") as f:
        lines = f.readlines()
    conv_pso = False
    sett=strip_setting_name(from_log_to_sett(log))
    for l in lines:
        if re.search("Converged",l):
            conv_pso = True
            conv_txt = f"\nPSO Converged after {str(l.split('after')[1])}"
            if verbose:
                print(conv_txt)
    if not conv_pso:
        sett=strip_setting_name(from_log_to_sett(log))
        conv_txt = f"\nWARNING: PSO of {sett} not converged in {log}\n"
        print(conv_txt)
    if save:
        with open(create_path_from_list([get_savefigpath(sett),"pso_convergence.txt"]),"a+") as f:
            f.write(f"{log}: {conv_txt}")

def check_pso_conv(setting,save=True,verbose=False):
    log = find_latest_log(setting)
    if verbose:
        print(f"Considering log file {log}")
    _check_pso_conv(log,save,verbose)

def check_success(settings,commands=[],verbose=0):   
    if type(settings)==str:
        settings=[settings]
    if commands==[]:
        from last_commands import progs as commands
    for i_set,setting_name in enumerate(settings):
        conv_pso = None    
        if "log_" ==setting_name[:4]:
            log = setting_name
            _check_pso_conv(log,save=True,verbose=True)
            setting_name=from_log_to_sett(log)
            if verbose>0:
                print("given log file ", log," considering setting file", setting_name)
        else:
            check_pso_conv(setting_name,save=True,verbose=True)
        if verbose>0:
            print(len(setting_name)*"#"+"######")
            print("# ",setting_name," #")
            print(len(setting_name)*"#"+"######\n")
        all_succ = 0
        for prog in commands:
            short_setting=strip_setting_name(setting_name)
            name_logs = "logs/"+short_setting+"_"+prog+"_*.log"
            list_of_logs = glob.glob(name_logs) 
            try:
                latest_log = max(list_of_logs, key=os.path.getctime)
            except ValueError:
                print("name_logs",name_logs)
                print("list_of_logs",list_of_logs)
                raise ValueError("Something is wrong, no log files found")
            if len(latest_log)==0:
                print(name_logs)
            success = False
            with open(latest_log,"r") as f:
                lines = f.readlines()
                last_lines = lines[-5:]
            for l in last_lines:
                if re.search("SUCCESS",l):
                    if verbose>1:
                        print("#"*80)
                        print("{ll:<55} successful".format(ll=latest_log))
                        print("#"*80)
                        print("\n")
                    success=True
                    all_succ+=1
                    break
            if not success:
                print("#"*80)
                print("WARNING {prg:<20} NOT SUCCESSFUL!\nLog: {ll}".format(prg=prog,ll=latest_log))
                print("\n")
                for l in last_lines:
                    print(l[:-1])
                print("\n"+"#"*80)
                print("\n")
                
        if all_succ==len(commands):
            stgn = strip_setting_name(setting_name)
            print("{sn:<20} ALL SUCCESSFUL".format(sn=stgn))
            
if __name__=="__main__":
    parser = ArgumentParser(description="Check if logs are successful")
    parser.add_argument("-CMD", "--commands",dest="commands", default=[],nargs="+",
                        help="Which commands have you run")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False,
                        help="Verbosity")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

    present_program(sys.argv[0])

    args = parser.parse_args()
    settings = args.SETTING_FILES
    commands = args.commands
    verbose = args.verbose
    
    if verbose:
        verbose = 2
    else:
        verbose = 1
    check_success(settings,commands,verbose)
    
    success(sys.argv[0])    

    
