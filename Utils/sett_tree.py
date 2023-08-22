#!/usr/bin/env python
# coding: utf-8

# I want to create a "genealogic tree" of the setting files
import sys
import graphviz
from os import walk
from argparse import ArgumentParser as ap

from Utils.tools import *

def find_EVERY_setting_path(setting_name):
    directories  = next(walk("."), (None, None, []))[1]
    setting_dir  = [i for i in directories if "setting" in i]
    setting_path = []
    for std in setting_dir:
        for j in next(walk("./"+std), (None, None, []))[2] :
            if setting_name in j:
                setting_path.append(std+"/"+j)
    if setting_name in next(walk("."), (None, None, []))[2]:
            setting_path.append(".")
    if setting_path==[]:
        raise RuntimeError("Setting file not found:"+setting_name)
    return setting_path
"""
def get_setting_from_path(setting_path):
    path_setting = str(setting_path).split("/")
    path = "/".join(path_setting[:-1])
    setting = path_setting[-1]
    sys.path.append(path)
    setting_module = importlib.import_module(setting) 
    return setting_module.setting()

def get_setting_with_path(setting,path):
    sys.path.append(path)
    setting_module = importlib.import_module(setting) 
    return setting_module.setting()
"""
def split_comment(comm):
    if ":" in comm:
        orig   = comm.split(":")[0]
        comm2  = " ".join(comm.split(":")[1:])
        orig_l = orig.split(" ")
        for i in orig_l:
            if "setting" in i or ".py" in i or any([filt in i for filt in ["160","475","814"]]):
                return i,comm2
    else:
        orig = comm.split(" ")
        for i in orig:
            if "setting" in i:
                return i,comm
    raise RuntimeError("No original setting name found:",comm)


def censure_comment(comm):
    # to avoid problems with HTML and special characters
    comm = comm.replace("&","_and_")
    comm = comm.replace("==","_ie_")
    comm = comm.replace("->","_to_")
    comm = comm.replace("=>","_to_")
    comm = comm.replace("''","_quote_")
    symbols = ("<",">","/","'","=")
    names   = ("_less than_","_larger than_","_slash_","_quote_","_equal_")
    for s,n in zip(symbols,names):
        comm = comm.replace(s,n)
    return comm

def redo_comment(comm):
    comm = censure_comment(comm)
    if len(comm)>80:
        _comm = comm.split(" ")
        _comm = [*_comm[:int(len(_comm)/2)],"<BR />",*_comm[int(len(_comm)/2):]]
        comm = " ".join(_comm)
    return comm
    
def get_string_node(name,comm):
    return "<<FONT POINT-SIZE='20'>"+name+ "</FONT> <BR /> <FONT POINT-SIZE='15'>"+comm+"</FONT>>"
    
    
if __name__=="__main__":
    #############################
    present_program(sys.argv[0])
    #############################
    parser = ap(description="Plot setting tree")
    #parser.add_argument("SETTING_FILES",nargs="+",default=[],help="setting file(s) to consider")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default="",
                    help="Directory name where to save the plot")
    parser.add_argument("-FM","--fig_mask", dest="fig_mask", default=False,action="store_true",
                        help="Add the masked input at the node")
    parser.add_argument("-cd","--common_denominator",type=str,dest="common_den", default="",
                    help="Common denominator (eg. filer: f105w, model: CP) for the naming of each setting file")
    args = parser.parse_args()
    dir_name      = args.dir_name
    common_den    = args.common_den
    fig_mask      = args.fig_mask
    #setting_names = args.SETTING_FILES  
    setting_names  = []
    backup_path = "backup_results/sett_tree_genealogy/"
    ##########################################################################################
    if common_den!="":
        if setting_names!=[]:
            raise RuntimeWarning("If common denominator is given, no setting files have to be considered")
        setts_paths = find_EVERY_setting_path(common_den)
        # the path should be only 1 (at least for now)
        test_path = [] 
        for paths in  setts_paths:
            test_path.append(paths.split("/")[0])
        #if len(list(set(test_path)))>1:
        #    raise RuntimeError("More then 1 directory for the settings is found! for now we expect only 1")
        path = test_path[0]
        setting_names = []
        for s in setts_paths:
                setting_names.append(verify_setting_name(s.split("/")[-1]))
    ##########################################################################################
    if dir_name=="":
        if common_den=="":
            for sett in setting_names:
                dir_name+=get_setting_name(sett)+"_"
            dir_name=dir_name[:-1]
        else:
            dir_name = common_den
    savedir = backup_path+"/"+dir_name
    mkdir(savedir)
    settings = []
    for sett in setting_names:
        try:
            settings.append(get_setting_module(sett,1))
        except:
            print("Setting ",sett," not up-to-date. Ignored") 
    s = graphviz.Digraph('settings-geanology')
    image = ""
    for sett in settings:
        comment   = sett.comments
        orig,comm = split_comment(comment)
        orig_name = strip_setting_name(orig)
        sett_name = strip_setting_name(sett)
        try:
            try:
                _,orig_comm = split_comment(get_setting_module(orig_name,1).comments)
            except RuntimeError:    
                print(orig_name," doesn't exist")
                orig_comm = "Inexistent"
        except AttributeError:
            print(orig_name," outdated")
            orig_comm = "Outdated"
        comm      = redo_comment(comm)
        orig_comm = redo_comment(orig_comm)
        str_sett  = get_string_node(sett_name,comm)
        str_orig  = get_string_node(orig_name,orig_comm)
        if fig_mask:
            savepath = get_savefigpath(sett_name).replace("./",os.getcwd()+"/")
            image = f"{savepath}/image_mask.png"
            fontcolor  = "white"
            if not os.path.exists(image):
                image=""
                fontcolor= "black"
        s.node(str_orig,label=str_orig)# "<<FONT POINT-SIZE='20'>"+orig_name+ "</FONT> >")
        s.node(str_sett,label=str_sett,image=image,imagescale="width",imagepos="mr",fontcolor=fontcolor)
        s.edge(str_orig,str_sett)
    #s = s.unflatten(stagger=6)  
    s.render(directory=savedir)
    print("Saved genealogy tree in "+savedir)
    success(sys.argv[0])

