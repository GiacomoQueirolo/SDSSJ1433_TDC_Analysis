# Basic functions 

import sys,os
import importlib
import json,pickle
from os import walk
import pathlib as pth
from datetime import datetime
from numpy import round as npround

def pycommand_info(date=True,toprint=False):
    # print the command and the time when the script was called 
    t_string = ""
    if date:
        t_string = " # "+str_date()
    cmd = "python "+" ".join(sys.argv)+t_string
    if toprint:
        print(cmd)
    return cmd 
    
def save_log_command(save_dir=".",date=True):
    # save the command line in a log file to relunch the program if needed
    str_date=""
    if date:
        if type(date) is str:
            str_date="_"+date
        else:
            str_date="_"+datetime.now().strftime("%d%m%Y")
    with open(save_dir+'/log_cmd'+str_date+'.txt', 'w') as f:
        f.write(pycommand_info())

def present_program(name_prog=""):
    if "./" == name_prog[:2]:
        name_prog = name_prog[2:]
    name_prog = str(name_prog).replace(".py","")
    line02 = "#"*int(len(name_prog)+4)
    line1 = "# "+name_prog+" #"
    print(line02)
    print(line1)
    print(line02)
    
def success(name_prog=""):
    line = "#"*int(len(name_prog)) + 12*"#"
    print(line)
    print(str(name_prog)+": SUCCESS!")
    print(line)
    
def separator(n_ash=100,n_lines=10):
    line = "#"*n_ash
    for i in range(n_lines):
        print(line)
        
def ProgressBar(i,max_,postText=""):
    n_bar =50 #size of progress bar
    j= i/max_
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%  {postText}")
    sys.stdout.flush()
"""
def is_defined_in_settings(thing,setting_name):
    # check if something is defined in the setting file
    # faster and more correct then "appears_in_setting"
    if type(thing) is not str:
        raise RuntimeWarning("Thing must be given as string")
    import importlib
    import sys
    path = find_setting_path(setting_name)
    if path!=".":
        sys.path.append(path)
    mod = importlib.import_module(setting_name)
    # Determine a list of names to copy to the current name space
    try:
        getattr(mod,thing)
        thing_is_defined = True
    except AttributeError:
        thing_is_defined = False
    return thing_is_defined
"""


def mkdir(dir_path,parents=True):
    pth.Path(dir_path).mkdir(parents=parents, exist_ok=True)
       
def say(something="lol"):
    import os
    os.system("say "+something)

def my_machine():
    from os import uname
    if uname().nodename=="jackquei-thinkpad":
        return True
    else:
        return False
        
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

        
def print_res_w_err(res,err,outstr=True):
    # print result with meaningfully rounded error
    exp=int(str("{:.1e}".format(err)).split("e")[-1])
    if exp>-2 and exp<1:
        exp=-1
    str_res = str(npround(res,abs(exp)))+" $\pm$ "+str(npround(err,abs(exp)))
    # to implement in case of very large/small numbers to use the scientific notation -> not the case here
    if outstr:
        return str_res
    else:
        print(str_res)
        
##########################
# setting file functions #
##########################

def _verify_setting_name(setting_name):
    # this way it works also with only the "f140w_SP_I" for example
    if type(setting_name) is str:
        setting_name = "settings_"+setting_name.replace("settings_","").replace(".py","")+".py"
        return setting_name
    else: 
        try:
            try: #if module
                setting_name.setting()
            except AttributeError: # if setting
                setting_name.image_name
        except: # pragma: no cover
            raise RuntimeError(str(setting_name)+" must be either str, module or setting, not "+str(type(setting_name)))
        return setting_name
        
def verify_setting_name(setting_name):
    # this way it works also with lists of settings
    if type(setting_name) is list:
        st_list= []
        for setting_name_i in setting_name: 
            st_list.append(_verify_setting_name(setting_name_i))
        return st_list
    else: 
        return _verify_setting_name(setting_name)
    
def find_setting_path(setting_name):
    setting_name =  verify_setting_name(setting_name)
    directories = next(walk("."), (None, None, []))[1]
    setting_dir = [i for i in directories if "setting" in i]
    setting_path = []
    for std in setting_dir:
        if setting_name in next(walk("./"+std), (None, None, []))[2]:
            setting_path.append(std)
    if setting_name in next(walk("."), (None, None, []))[2]:
            setting_path.append(".")
    if setting_path==[]:
        raise RuntimeError("Setting file not found:"+setting_name)
    elif len(setting_path)>1:
        string= "More then one setting found with same name in:"
        for st in setting_path:
            string+=st+"\n"
        raise RuntimeError(string)
    else:
        return setting_path[0]

def check_if_CP(setting):
    if type(setting) is str:
        setting=get_setting_module(setting).setting()
    try:
        CP = setting.CP
    except AttributeError:
        if "gamma" in list(setting.lens_params[0][0].keys()):
            CP=True
        else:
            CP=False
    return CP

def check_if_WS(setting):
    if type(setting) is str:
        setting=get_setting_module(setting).setting()
    try:
        WS = setting.WS
    except AttributeError:
        if hasattr(setting,"source_params"):
                WS=False
        else:
                WS=True
    return WS

def check_if_SUB(setting):
    setting = get_setting_module(setting).setting()
    if setting.sub is True:
        # no matter the lens_light_model_name: it could be None bc it was already subtracted
        SUB = True
    else:
        SUB = False
    #else: 
    #    raise RuntimeError("Something wrong:"+str(setting.lens_light_model_name)+" "+str(setting.sub))
    return SUB


def _get_setting_module(setting_name,sett=False):
    setting_name = verify_setting_name(setting_name)
    try:
        setting_name.setting()
        if sett:        
            return setting_name.setting()
        else:
            return setting_name
    except AttributeError:
        try:
            setting_name.comments
            if sett:
                return setting_name
            else:
                raise AttributeError
        except AttributeError:     
            setting_name     = get_setting_name(setting_name)
            setting_position = find_setting_path(setting_name)
            sys.path.append(setting_position)
            setting_module   = importlib.import_module(setting_name.replace(".py","")) 
            if sett:
                return setting_module.setting()
            else:
                return setting_module

def get_setting_module(setting_name,sett=False):
    if type(setting_name) is list:
        stmd = []
        for s in setting_name:
            stmd.append(_get_setting_module(s,sett=sett))
        return stmd
    else:
        return _get_setting_module(setting_name,sett=sett)
    
def get_savefigpath(setting,backup_path="backup_results"):
    backup_path  = str(backup_path)    
    setting_name = get_setting_name(setting)
    setting_name = setting_name.replace(".py","") 
    savefig_path = "./"+backup_path+"/"+setting_name.replace("settings_","")+"/"
    return savefig_path

def get_savemcmcpath(setting,backup_path="backup_results"):
    backup_path   = str(backup_path)    
    setting_name  = get_setting_name(setting)
    setting_name  = setting_name.replace(".py","")  
    savemcmc_path = "./"+backup_path+"/"+setting_name.replace("settings_","mcmc_")+"/"
    return savemcmc_path

def get_backend_filename(setting,backup_path="backup_results"):
    setting_name  = get_setting_name(setting)
    savemcmc_path = get_savemcmcpath(setting,backup_path=backup_path)
    backend_name  = create_path_from_list([savemcmc_path,
                    "backend_mcmc_"+strip_setting_name(setting_name)+".hp5"])
    return backend_name
    

    
def get_kwargs_params(setting_name):
    setting_mod = get_setting_module(setting_name)
    sett = setting_mod.setting()
    kwargs_params = {'lens_model': sett.lens_params,
                'point_source_model': sett.ps_params,
                'lens_light_model': sett.lens_light_params}
    try:
        kwargs_params['source_model'] = sett.source_params
    except AttributeError:
        # it's a _ws setting file
        pass
    return kwargs_params
    
def print_kwargs(setting_name):
    params_type = ["init","sigma","fixed","min","max"]
    kw = get_kwargs_params(setting_name)
    print("")
    for key in list(kw.keys()):
        model = kw[key]
        print(len(key)*"#"+"####")
        print("#",key,"#")
        print(len(key)*"#"+"####\n")
        for i_m,mod in enumerate(model):
            print(params_type[i_m])
            print(len(params_type[i_m])*"#","\n")
            for prof in mod:
                print(prof)
            print("")
        print("")

def _get_setting_name(setting):
    setting = verify_setting_name(setting)
    if type(setting) is not str:
        try:
            setting_name = setting.__module__
        except AttributeError:
            setting_name = setting.setting().__module__
    else:
        setting_name = setting
    return str(setting_name)
    
def get_setting_name(setting):
    if type(setting) is list:
        st_nm = []
        for sett_i in setting:
            st_nm.append(_get_setting_name(sett_i))
        return st_nm
    else:
        return _get_setting_name(setting) 
    
def strip_setting_name(setting,filter=False):
    setting_name = get_setting_name(setting)
    strip        = setting_name.replace("settings_","").replace(".py","")
    if filter:
        try:
            strip = get_setting_module(setting,1).filter_name
        except:
            strip = strip.split("_")[0]
    return strip

def print_setting(setting):
    sett = strip_setting_name(setting)
    print(len(sett)*"#"+"####")
    print(sett)
    print(len(sett)*"#"+"####")

def get_filter(setting):
    return strip_setting_name(setting,filter=True)

def create_path_from_list(ordered_list):
    # we asume an ordered list of directories and subdirectories (and maybe a file in the end)
    path=""
    for p in ordered_list:
        if type(p) is not str:
            raise RuntimeError(f"Input must be a list of str, not {type(p)}")
        if p[0]!="/" and p[0]!=".":
            p=f"/{p}"
        path+=p
    return path
        
def save_json(data,filename):
    from numpy import array 
    data = array(data).tolist()
    with open(filename, 'w') as f: 
        json.dump(data, f)

def save_json_name(setting,path,filename):
    setting_name = get_setting_name(setting)
    json_name = setting_name.replace("settings",filename).replace(".py","")+".json"
    name = create_path_from_list([path,json_name]) 
    return name


def save_mcmc_json(setting,data,filename,backup_path="backup_results"):
    setting_name  = get_setting_name(setting)
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
    name = save_json_name(setting,savemcmc_path,filename)
    save_json(data,name)
    
def create_dir_name(settings,save_dir=".",dir_name=None,backup_path="./backup_results",copy_settings=True):
    save_dir = str(backup_path)+"/"+str(save_dir)+"/"
    if type(settings) is not list:
        settings=[settings]
    if dir_name is None:
        filters=[strip_setting_name(st) for st in settings]
        lf = ""
        for f in filters:
            lf+=f+"_"
        lf = lf[:-1]
        save_dir=save_dir+"/"+str(lf)
    else:
        save_dir=save_dir+"/"+dir_name
    mkdir(save_dir)
    if copy_settings:
        setting_position = []
        for sn in settings:
            setting_position.append(find_setting_path(sn))
        for st_i in range(len(setting_position)):
            os.system("cp "+str(setting_position[st_i])+"/"+str(get_setting_name(settings[st_i]))+    " "+str(save_dir)+"/.")

    return save_dir


def get_param_input(settings,param_name,prm_type=None):
    sett = get_setting_module(settings).setting()
    if "image" in param_name:
        prms = sett.ps_params
        try:
            int(param_name[-1])
        except ValueError:
            raise ValueError("Give also the number of the image")
            
    if "lens" in param_name:
        if "light" in param_name:
            prms = sett.lens_light_params
        else:
            prms = sett.lens_params
    if "source" in param_name:
        if check_if_WS(sett):
            print("WARNING: This setting file doesn't have a source, hence no ",param_name)
            return -1
        prms = sett.source_params
    
    prm_nm  = "_".join(param_name.split("_")[:-1]) # get rid of the type of param 
    prm_nmb = int(param_name[-1])
    kwg_prms = {}
    if not "image" in param_name:
        try:
            kwg_prms["fixed"] = prms[2][prm_nmb][prm_nm]
        except KeyError:
            kwg_prms["init"]  = prms[0][prm_nmb][prm_nm]
            kwg_prms["sigma"] = prms[1][prm_nmb][prm_nm]
            kwg_prms["lower"] = prms[3][prm_nmb][prm_nm]
            kwg_prms["upper"] = prms[4][prm_nmb][prm_nm]
    else: 
        try:
            kwg_prms["fixed"] = prms[2][0][prm_nm][prm_nmb]
        except KeyError:
            kwg_prms["init"]  = prms[0][0][prm_nm][prm_nmb]
            kwg_prms["sigma"] = prms[1][0][prm_nm][prm_nmb]
            kwg_prms["lower"] = prms[3][0][prm_nm][prm_nmb]
            kwg_prms["upper"] = prms[4][0][prm_nm][prm_nmb]
    if prm_type:
        return kwg_prms[prm_type]
    else:
        return kwg_prms

def setting_is_infrared(setting):
    filt = get_filter(setting).lower()
    if "f1" in filt:
        return True
    else:
        return False

def str_date():
    now      = datetime.now()
    t_string = now.strftime("%d/%m/%Y %H:%M:%S")
    return t_string

def pickle_results(res,name,savefig_path=""):
    # I save the kwargs result in a pickle, readable way
    if not ".data" in name:
        name+=".data"
    with open(savefig_path+name,"wb") as f:
        pickle.dump(res, f)

###
# for nice formatting:
import argparse
class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass 
###

# for pycs3
##########################################
def get_config(lensname,dataname=None,config_path="./myconfig/"):
    try:
        lensname.timeshifts
        return lensname # it is already the config file
    except AttributeError:
        if dataname is None:
            lensname,dataname=lensname.replace("myconfig_","").replace(".py","").split("_")
        sys.path.append(config_path+"/")
        config_file = "myconfig_" + lensname + "_" + dataname
        config = importlib.import_module(config_file)
        return config

def get_simdir(path):
    return str(path).replace("Analysis","Simulation").replace("analysis","simulation")

def get_analdir(path):
    return str(path).replace("Simulation","Analysis").replace("simulation","analysis")

# MOD_ROB
def check_success_analysis(path):
    #just to be sure:
    path    = get_analdir(path)
    success = pickle.load(open(path+"/success.data","rb"))
    return success["success"]
##########################################
