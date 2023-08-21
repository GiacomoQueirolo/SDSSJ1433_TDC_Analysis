import dill
from Utils.tools import *

def verify_multifilter_setting_name(multifilter_setting_name):
    # Setting_Multifilter is assumed that it can't be a list
    # this returns the name with which the multifilter is stored
    if type(multifilter_setting_name) is str:
        multifilter_setting_savename = f'mltf_setting_{multifilter_setting_name.replace("mltf_setting_","").replace(".dll","")}.dll'
        return multifilter_setting_savename
    else: 
        try:
            multifilter_setting_name.settings
        except: # pragma: no cover
            raise RuntimeError(f"{multifilter_setting_name} must be either str or Setting_Multifilter, not {type(multifilter_setting_name)}")
        return f'mltf_setting_{multifilter_setting_name.name}.dll'


def find_multifilter_setting_path(multifilter_setting,main_dir="./"):
    multifilter_setting_savename =  verify_multifilter_setting_name(multifilter_setting)
    directories = next(walk(main_dir), (None, None, []))[1]
    mltf_setting_dir = [i for i in directories if "multifilter" in i]
    mltf_setting_path = []
    for mstd in mltf_setting_dir:
        if multifilter_setting_savename in next(walk(main_dir+"/"+mstd), (None, None, []))[2]:
            mltf_setting_path.append(mstd)
    if multifilter_setting_savename in next(walk(main_dir), (None, None, []))[2]:
            mltf_setting_path.append(main_dir)
    if mltf_setting_path==[]:
        raise RuntimeError(f"Multifilter Setting file not found: {multifilter_setting_savename}")
    elif len(mltf_setting_path)>1:
        string= "More then one setting found with same name in:"
        for mst in mltf_setting_path:
            string+=mst+"\n"
        raise RuntimeError(string)
    else:
        return mltf_setting_path[0]

def get_multifilter_setting_name(multifilter_setting):
    # return the pure name of the multifilter setting file
    # !: different form get_setting_name, as that returns "settings_.."
    try:
        multifilter_setting.name
        return multifilter_setting.name
    except AttributeError:
        multifilter_setting_savename = verify_multifilter_setting_name(multifilter_setting)
        multifilter_setting_name     = multifilter_setting_savename.replace("mltf_setting_","").replace(".dll","")
        return str(multifilter_setting_name)
    
def _string_is_filter(str_i):
    str_i = str(str_i).replace("_","").lower()
    if len(str_i)!=5:
        return False
    if str_i[0]=="f":
        if str_i[-1]=="w" or str_i[-1]=="x":
            try:
                int(str_i[1:-1])
                return True
            except ValueError:
                return False
        
def get_settings_from_multifilter(multifilter_setting):
    multifilter_setting_savename  = verify_multifilter_setting_name(multifilter_setting).replace("mltf_setting_","").replace(".dll","")
    try:
        multifilter_setting.name
        return multifilter_setting.settings
    except AttributeError:
        split_name = multifilter_setting_savename.split("_")
        filt_name  = [_string_is_filter(s) for s in split_name]
        # discard what is present before the 1st filter
        split_name = split_name[filt_name.index(True):]
        filt_name  = [_string_is_filter(s) for s in split_name]
        sett_list   = []
        for i,(fn, sn) in enumerate(zip(filt_name,split_name)):
            if fn==True:
                try:
                    sett_list.append(st)
                except:
                    n=0
                st=[]
            st.append(sn)
            if i== len(filt_name)-1:
                sett_list.append(st)
        settings = []
        for st in sett_list:
            settings.append(get_setting_module("_".join(st),1))
        return settings
        
def get_multifilter_setting_module(multifilter_setting,main_dir="./"):
    multifilter_setting_savename  = verify_multifilter_setting_name(multifilter_setting)
    multifilter_setting_position  = find_multifilter_setting_path(multifilter_setting_savename,main_dir=main_dir)
    while True:
        try:
            with open(f"{multifilter_setting_position}/{multifilter_setting_savename}","rb") as f:
                multifilter_setting_module = dill.load(f)
            break
        except ModuleNotFoundError:
            nothing=True
        # we have to load the relative setting files first:
        for setting in get_settings_from_multifilter(multifilter_setting):
            get_setting_module(setting,1)
    return multifilter_setting_module
    
def check_mltf_setting(funct):
    def _check_mltf_sett(mltf_setting,*args,**kwargs):
        # note: mltf_setting MUST BE the first argument
        mltf_setting=get_multifilter_setting_module(mltf_setting)
        return funct(mltf_setting,*args,**kwargs)
    return _check_mltf_sett


def _append(kwargs,sett_kwargs):
    for sk in sett_kwargs:
        kwargs.append(sk)
    return kwargs

def strip_multifilter_name(multifilter_setting,filter=False):
    mltf_setting_name = get_multifilter_setting_name(multifilter_setting)
    strip             = mltf_setting_name.replace("mltf_setting_","").replace(".dll","")
    if filter:
        try:
            strip = "_".join(get_multifilter_setting_module(mltf_setting_name).filters)
        except:
            # pragma no cover
            raise RuntimeError("no filter found")
    return strip


# implemented in the Setting_multifilter class 
"""
# change of idea: in the directory they will be called only by their name, w/o
# reference to the name of the mltf (which is too long)
#def save_json_name_mltf(multifilter_setting,path,filename):
    # mlf_name   = get_multifilter_setting_name(multifilter_setting)
    #json_name  = mlf_name.replace("mltf_setting",filename).replace(".dill","")+".json"
def save_json_name_mltf(filename,path):
    json_name  = filename.replace(".dill","")+".json"
    name       = create_path_from_list([path,json_name]) 
    return name

"""