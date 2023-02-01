from argparse import ArgumentParser
import re
import glob
import os
from tools import present_program
import sys

parser = ArgumentParser(description="Check if logs are successful")
parser.add_argument("-CMD", "--comands",dest="comands", default=[],nargs="+",
                    help="Which comands have you run")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False,
                    help="Verbosity")
parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

present_program(sys.argv[0])

args = parser.parse_args()
settings = args.SETTING_FILES
comands = args.comands
verbose = args.verbose
for i_set,setting_name in enumerate(settings):    
    if "log_" ==setting_name[:4]:
        from_log_to_sett = setting_name[4:]
        if ".dat" == from_log_to_sett[-4:]:
            from_log_to_sett=from_log_to_sett[:-4]
        #get rid of date
        from_log_to_sett = "_".join(from_log_to_sett.split("_")[:-1])
        from_log_to_sett = "settings_"+from_log_to_sett+".py"
        settings[i_set] = from_log_to_sett
        print("given log file ", setting_name," considering setting file", from_log_to_sett)

if comands==[]:
    from last_comands import progs as comands
for setting_name in settings:
    print(len(setting_name)*"#"+"######")
    print("# ",setting_name," #")
    print(len(setting_name)*"#"+"######\n")
    all_succ = 0
    for prog in comands:
        short_setting=setting_name.replace("settings_","").replace(".py","")
        name_logs = "logs/"+short_setting+"_"+prog+"_*.log"
        list_of_logs = glob.glob(name_logs) 
        latest_log = max(list_of_logs, key=os.path.getctime)
        success = False
        with open(latest_log,"r") as f:
            lines = f.readlines()
            last_lines = lines[-5:]
        for l in last_lines:
            if re.search("SUCCESS",l):
                if verbose:
                    print(latest_log," successful")
                    print("\n")
                success=True
                all_succ+=1
                break
        if not success:
            print("WARNING ",prog," NOT SUCCESSFUL!\nLog: "+latest_log)
            print("\n")
            for l in last_lines:
                print(l[:-1])
    
    if all_succ==len(comands):
        print(setting_name, " ALL SUCCESSFUL")
