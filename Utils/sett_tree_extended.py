#!/usr/bin/env python
# coding: utf-8

# In[1]:


# I want to create a "genealogic tree" of the setting files

import sys
import graphviz
from os import walk
from PIL import Image
from argparse import ArgumentParser as ap


from Utils.tools import *
from Utils.sett_tree import *

def _resize_image(image_path,image_name,end_path,new_image_name=None,thumbnail=True):
    image = Image.open(image_path+"/"+image_name)
    px = np.asarray(image.convert("L"))
    coord_k = np.column_stack(np.where(px<10))
    edge0,edge1 = coord_k[0],coord_k[-1] 
    dx,dy = np.float(edge1[0]-edge0[0]),np.float(edge1[1]-edge0[1]) 
    shift = np.array([int(dx/4.),int(dy/4.)]) 
    px0,px1= np.array(edge0)+shift,np.array(edge1)-shift 
    shift_right = 20
    image_cropped = image.crop(box=(px0[0]+shift_right,px0[1],px1[0]+shift_right,px1[1]))
    if thumbnail:
        image.thumbnail((200,200),Image.ANTIALIAS)
    if not new_image_name:
        image_cropped.save(end_path+"/"+image_name)
    else:
        image_cropped.save(end_path+"/"+new_image_name)
        
if __name__=="__main__":
    #############################
    present_program(sys.argv[0])
    #############################
    parser = ap(description="Plot setting tree",formatter_class=CustomFormatter)
    #parser.add_argument("SETTING_FILES",nargs="+",default=[],help="setting file(s) to consider")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default="",
                    help="Directory name where to save the plot. If not given, same as Common denominator")
    parser.add_argument("-FM","--fig_mask", dest="fig_mask", default=False,action="store_true",
                        help="Add the masked input at the node")
    parser.add_argument("-cd","--common_denominator",type=str,dest="common_den",default=None,
                    help="Common denominator (eg. filer: f105w, model: CP) for the naming of each setting file")
    args = parser.parse_args()
    dir_name      = args.dir_name
    common_den    = args.common_den
    if common_den is None:
        raise RuntimeError("Give at least a common denominator")
    fig_mask      = args.fig_mask
    #setting_names = args.SETTING_FILES  
    setting_names  = []
    backup_path = "backup_results/sett_tree_genealogy/"
    ##########################################################################################
    if common_den!="":
        setts_paths = find_EVERY_setting_path(common_den)
        # the path should be only 1 (at least for now)
        test_path = [] 
        for paths in  setts_paths:
            test_path.append(paths.split("/")[0]) 
        path = test_path[0]
        setting_names = []
        for s in setts_paths:
                setting_names.append(verify_setting_name(s.split("/")[-1]))
    ##########################################################################################
    if dir_name=="":
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
    image_n = 0
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
        comm      = redo_comment(comm,25)
        orig_comm = redo_comment(orig_comm,25)
        str_sett  = get_string_node(sett_name,comm)
        str_orig  = get_string_node(orig_name,orig_comm)
        s.node(str_orig,label=str_orig)
        fontcolor= "black"
        if fig_mask:
            savepath = get_savefigpath(sett_name).replace("./",os.getcwd()+"/")
            image = f"{savepath}/image_mask.png"
            if not os.path.exists(image):
                s.node(str_sett,label=str_sett,fontcolor=fontcolor)
            else:
                fontcolor  = "white"
                new_image_name="image_"+str(image_n)+".png"
                image_n+=1
                _resize_image(savepath,"image_mask.png",savedir,new_image_name=new_image_name)
                image = f"{os.getcwd()}/{savedir}/{new_image_name}"
                s.node(str_sett,label=str_sett,image=image,imagepos="mc",fontcolor=fontcolor)
        else:
            s.node(str_sett,label=str_sett,fontcolor=fontcolor)
        s.edge(str_orig,str_sett) 
    s.render(directory=savedir)
    print("Saved genealogy tree in "+savedir)
    success(sys.argv[0])

