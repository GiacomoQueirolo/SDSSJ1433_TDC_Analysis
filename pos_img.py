#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# produce .reg file with the initial position of images and their freedom
import numpy as np
import sys
import argparse

from tools import *
from conversion import conv_radec_to_xy,qphi_from_e1e2
from image_manipulation import get_numPix,load_fits, get_rotangle
from source_pos import *
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


# In[ ]:


def radius(list_max,list_min,sets):
    max_ra,max_dec = list_max
    min_ra,min_dec = list_min
    max_x,max_y    = conv_radec_to_xy(sets,ra=max_ra,dec=max_dec)
    min_x,min_y    = conv_radec_to_xy(sets,ra=min_ra,dec=min_dec)
    dx,dy          = np.abs(np.array(max_x)-np.array(min_x)),np.abs(np.array(max_y)-np.array(min_y))
    
    rd_i       = [np.sqrt(dx_i**2+dy_i**2)/(2.) for dx_i,dy_i in zip(dx,dy)]
    if len(rd_i)==1:
        return rd_i[0]
    return rd_i

def get_dxdy(list_max,list_min,sets):
    max_ra,max_dec = list_max
    min_ra,min_dec = list_min
    dra,ddec       = np.abs(np.array(max_ra)-np.array(min_ra)),np.abs(np.array(max_dec)-np.array(min_dec))
    return(dra/sets.pix_scale,ddec/sets.pix_scale)

"""
# WIP
def get_boundaries(cnt,points):
    # given a center and a set of 2D points
    # define the boundaries by selecting the "most extreme points" 
    # for each direction, ie direction between the center and the points
    xc,yc = cnt
    angles = []
    boundary = []
    for p in points:
        xp,yp = p
        dx,dy = xp-xc,yp-yc
        d = np.sqrt(dx**2+dy**2)
        angles.append(np.arctan2(dy,dx))
    
    angle_sorted = sorted(angles)
    diff_angle_sorted = [a_s[ias]-a_s[ias-1] for ias in range(len(angle_sorted[1:]))]
    
    thresh = 0.001
    i = 0
    boundary = [points]
    while i<len(points):
        diff_angle_sorted = 
        
        
    for das in diff_angle_sorted:
        if das < thresh:
            
        
        add = True
        for i_a,a in enumerate(angles):
            if abs(a-angle_i)<0.001:
                same_angle=True
                other_p = points[i_a]
                othxp,othyp = other_p
                othdx,othdy = othxp-xc,othyp-yc
                other_d = np.sqrt(othdx**2+othdy**2)
                if d>other_d:
                    add=True
                else:
                    add=False
                    break
        angles.append(angle_i)
        if add:
            boundary.append(p)
    return boundary
"""
    
def create_pos_img(setting_name,RP,RM,also_fits=False,save=False):
    sets    = get_setting_module(setting_name).setting()
    numPix  = get_numPix(sets)
    
    #init images
    maxra,maxdec = sets.ps_params[-1][0]["ra_image"],sets.ps_params[-1][0]["dec_image"]
    minra,mindec = sets.ps_params[-2][0]["ra_image"],sets.ps_params[-2][0]["dec_image"]
    rd_i         = radius([maxra,maxdec],[minra,mindec],sets)
    dx_i,dy_i    = get_dxdy([maxra,maxdec],[minra,mindec],sets)
    centra,centdec = sets.ps_params[0][0]["ra_image"],sets.ps_params[0][0]["dec_image"]
    centx,centy    = conv_radec_to_xy(sets,ra=centra,dec=centdec)
    
   
    # write images
    let = ["A","C","B","D"]
    # to_plot images
    im_lbl = ["Init. image pos. (Rad approx. Max Freedom)",None,None,None]
    to_plot =  [{"x":centx[i],"y":centy[i], "rds": rd_i[i], "txt": let[i],"col":"lime","lns":"-", "lbl":im_lbl[i]} for i in range(len(centx))]
    
    if also_fits:
        str_reg = '# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nimage\n'
        coord_y = 4
        for i in range(len(centx)):
            #str_reg+="circle("+str(centx[i])+",	"+str(centy[i])+", "+str(rd_i[i])+") # color=green  text={"+str(let[i])+"} \n"    
            str_reg+="box("+str(centx[i]+1)+",	"+str(centy[i]+1)+", "+str(dx_i[i])+","+str(dy_i[i])+",0) # color=green  text={"+str(let[i])+"} \n"    
        str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Initial image pos. (Rad=Max Freedom) } color=lime\n"
        coord_y +=2
        
    # init source #
    kw_lens = sets.lens_params[0]
    source_ra_max_,source_dec_max_ = get_source_gen(maxra,maxdec,kw_lens,sets)
    source_ra_min_,source_dec_min_ = get_source_gen(minra,mindec,kw_lens,sets)
    source_ra_max,source_dec_max   = np.max(source_ra_max_),np.max(source_dec_max_)
    source_ra_min,source_dec_min   = np.min(source_ra_min_),np.min(source_dec_min_)
    rd_src                         = radius([source_ra_max,source_dec_max],[source_ra_min,source_dec_min],sets)
    dx_src,dy_src                  = get_dxdy([source_ra_max,source_dec_max],[source_ra_min,source_dec_min],sets)
    source_ra_cnt_,source_dec_cnt_ = get_source_gen(centra,centdec,kw_lens,sets)
    source_ra_cnt,source_dec_cnt   = np.no more needed to correct for DS9 starting from 1mean(source_ra_cnt_),np.mean(source_dec_cnt_)
    source_x_cnt,source_y_cnt      = conv_radec_to_xy(setting_name,ra=source_ra_cnt,dec=source_dec_cnt)
    if also_fits:
        ## add source ##
        str_reg+="box("+str(source_x_cnt[0]+1)+", "+str(source_y_cnt[0]+1)+", "+str(dx_src)+" , "+str(dy_src)+",0. ) # color=yellow  text={Source} \n"
        #str_reg+="circle("+str(source_x_cnt[0])+",	"+str(source_y_cnt[0])+", "+str(rd_src)+") # color=yellow  text={Source} \n"
        str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Initial source pos. (Rad=Max Freedom) } color=yellow \n"
        coord_y+=2    
    # to_plot src
    to_plot.append({"x":source_x_cnt[0],"y":source_y_cnt[0], "rds":rd_src, "txt":"Source","col":"y","lns":"-", "lbl":"Init. source pos. (Rad approx. Max Freedom)"})

    # main lens
    maxra_lns,maxdec_lns = sets.lens_params[-1][0]["center_x"],sets.lens_params[-1][0]["center_y"]
    minra_lns,mindec_lns = sets.lens_params[-2][0]["center_x"],sets.lens_params[-2][0]["center_y"]
    rd_lns               = radius([maxra_lns,maxdec_lns],[minra_lns,mindec_lns],sets)
    dx_lns,dy_lns        = get_dxdy([maxra_lns,maxdec_lns],[minra_lns,mindec_lns],sets)
    centra_lns,centdec_lns = sets.lens_params[0][0]["center_x"],sets.lens_params[0][0]["center_y"]
    centx_lns,centy_lns    = conv_radec_to_xy(setting_name,ra=centra_lns,dec=centdec_lns)
    
    ## add main_lens ##
    if also_fits:
        #str_reg+="circle("+str(centx_lns[0])+",	"+str(centy_lns[0])+", "+str(rd_lns)+") # color=white  text={Main lens} \n"
        str_reg+="box("+str(centx_lns[0]+1)+",	"+str(centy_lns[0]+1)+", "+str(dx_lns)+", "+str(dy_lns)+", 0.) # color=white  text={Main lens} \n"
        str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Initial lens pos. (Rad=Max Freedom) } color=white \n"
        coord_y+=2
    # to_plot main lens
    to_plot.append({"x":centx_lns[0],"y":centy_lns[0], "rds":rd_lns, "txt":"Main Lens","col":"w","lns":"-", "lbl":"Init. main lens pos. (Rad approx. Max Freedom)"})

    # perturber
    maxra_prt,maxdec_prt = sets.lens_params[-1][1]["center_x"],sets.lens_params[-1][1]["center_y"]
    minra_prt,mindec_prt = sets.lens_params[-2][1]["center_x"],sets.lens_params[-2][1]["center_y"]
    rd_prt               = radius([maxra_prt,maxdec_prt],[minra_prt,mindec_prt],sets)
    dx_prt,dy_prt        = get_dxdy([maxra_prt,maxdec_prt],[minra_prt,mindec_prt],sets)
    centra_prt,centdec_prt = sets.lens_params[0][1]["center_x"],sets.lens_params[0][1]["center_y"]
    centx_prt,centy_prt    = conv_radec_to_xy(setting_name,ra=centra_prt,dec=centdec_prt)
    ## add perturber ##
    if also_fits:
        #str_reg+="circle("+str(centx_prt[0])+",	"+str(centy_prt[0])+", "+str(rd_prt)+") # color=violet  text={Pert. lens} \n"
        str_reg+="box("+str(centx_prt[0]+1)+", "+str(centy_prt[0]+1)+", "+str(dx_prt)+" , "+str(dy_prt)+",0. ) # color=violet  text={Pert. lens} \n"
        str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Initial perturber pos. (Rad=Max Freedom) } color=violet \n"
        coord_y+=2
    # to_plot main lens
    to_plot.append({"x":centx_prt[0],"y":centy_prt[0], "rds":rd_prt, "txt":"Perturber","col":"violet","lns":"-", "lbl":"Init. perturber pos. (Rad approx. Max Freedom)"})
    
        
    if RP:
        from source_pos import get_source_pos_PSO
        kw_res = get_kwres(setting_name,updated=False)["kwargs_results"]
        ps_res = kw_res["kwargs_ps"][0]
        ra_res,dec_res = ps_res["ra_image"],ps_res["dec_image"]
        X_res,Y_res    = conv_radec_to_xy(setting_name,ra=ra_res,dec=dec_res)
        if also_fits:
            for i in range(len(X_res)):
                #str_reg+="circle("+str(X_res[i])+",	"+str(Y_res[i])+", "+str(rd_i[i]/2.)+") # color=white dash=1\n"
                str_reg+="point("+str(X_res[i]+1)+",	"+str(Y_res[i]+1)+")# point=x 20 color=cyan dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Image pos. from PSO } color=cyan\n"
            coord_y+=2
        for i in range(len(X_res)):
            rd = to_plot[i]["rds"]/2.
            if i==0:
                lbl_im = "Image pos. from PSO"
            else:
                lbl_im = None
            to_plot.append({"x":X_res[i] ,"y":Y_res[i], "rds":rd,"txt":None,"col":"cyan","lns":"cross", "lbl":lbl_im})
        
        #### add source ###
        kw_res_src     = get_source_pos_PSO(sets)
        ra_src,dec_src = kw_res_src["source_ra"],kw_res_src["source_dec"]
        X_src,Y_src    = conv_radec_to_xy(setting_name,ra=ra_src,dec=dec_src)
        #str_reg+="circle("+str(X_src[0])+",	"+str(Y_src[0])+","+str(rd_src/2.)+") # color=cyan dash=1\n"
        if also_fits:
            str_reg+="point("+str(X_src[0]+1)+",	"+str(Y_src[0]+1)+") # point=x 20 color=orange dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Source pos. from PSO } color=orange\n"
            coord_y+=2        
        rd = to_plot[4]["rds"]/2.
        to_plot.append({"x":X_src[0],"y":Y_src[0], "rds":rd,"txt":None,"col":"orange","lns":"cross", "lbl":"Source pos. from PSO"})
        ##### add main lens #####
        kw_lens        = kw_res["kwargs_lens"]
        ra_lns,dec_lns = kw_lens[0]["center_x"],kw_lens[0]["center_y"]
        X_lns,Y_lns    = conv_radec_to_xy(setting_name,ra=ra_lns,dec=dec_lns)
        if also_fits:
            str_reg+="point("+str(X_lns[0]+1)+",	"+str(Y_lns[0]+1)+") # point=x 20 color=magenta dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Main lens pos. from PSO } color=magenta\n"
            coord_y+=2
        
        to_plot.append({"x":X_lns[0],"y":Y_lns[0], "rds":10,"txt":None,"col":"magenta","lns":"cross", "lbl":"Main lens pos. from PSO"})
        #### add vector indicating pointing angle of Lens Mass (only on fits) ####
        if also_fits:
            e1,e2   = kw_lens[0]["e1"],kw_lens[0]["e2"]
            _,phi   = qphi_from_e1e2(e1,e2,ret_deg=True)
            # get limits of phi
            e1_min,e1_max = sets.lens_params[-2][0]["e1"],sets.lens_params[-1][0]["e1"]
            e2_min,e2_max = sets.lens_params[-2][0]["e2"],sets.lens_params[-1][0]["e2"]
            e1_prior = np.random.uniform(e1_min,e1_max,50000)
            e2_prior = np.random.uniform(e2_min,e2_max,50000)
            _,phi_prior = qphi_from_e1e2(e1_prior,e2_prior,ret_deg=True)
            phi_max,phi_min = phi_prior.max(), phi_prior.min()
            print("phi_max,phi_min,phi",phi_max,phi_min,phi)
            print("e1_min,e1_max,e1",e1_min,e1_max,e1)
            print("e2_min,e2_max,e2",e2_min,e2_max,e2)
            plt.hist(np.array(phi_prior),histtype="step")
            plt.axvline(phi,label="measured phi from PSO")
            plt.legend()
            plt.savefig(get_savefigpath(sets)+"/ellipticity/phi_vs_prior.png")
            plt.close()
            rotang  = get_rotangle(sets) 
            phi_rt      = - phi + rotang 
            phi_rt_max  = - phi_max + rotang 
            phi_rt_min  = - phi_min + rotang 
            lenght  = 10 

            str_reg+="# vector("+str(X_lns[0]+1)+",	"+str(Y_lns[0]+1)+"," +str(lenght)+"," + str(phi_rt)+")  vector=1 color=springgreen dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Main lens P.A. from PSO } color=springgreen\n"
            coord_y+=2
            str_reg+="# vector("+str(X_lns[0]+1)+",	"+str(Y_lns[0]+1)+"," +str(lenght)+"," + str(phi_rt_max)+")  vector=1 color=tomato dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Max main lens P.A. from prior } color=tomato\n"
            coord_y+=2
            str_reg+="# vector("+str(X_lns[0]+1)+",	"+str(Y_lns[0]+1)+"," +str(lenght)+"," + str(phi_rt_min)+")  vector=1 color=cornflowerblue dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Min main lens P.A. from prior} color=cornflowerblue\n"
            coord_y+=2

            if not check_if_SUB(sets):
                kw_lens_light = kw_res["kwargs_lens_light"]
                e1_ll,e2_ll = kw_lens_light[0]["e1"],kw_lens_light[0]["e2"]
                _,phi_ll    = qphi_from_e1e2(e1_ll,e2_ll,ret_deg=True)
                phi_rt_ll   = - phi_ll + rotang
                str_reg+="# vector("+str(X_lns[0]+1)+",	"+str(Y_lns[0]+1)+"," +str(lenght)+"," + str(phi_rt_ll)+")  vector=1 color=gold \n"
                str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Main lens light P.A. from PSO } color=gold\n"
                coord_y+=2
        
        ##### add perturber #####
        ra_prt,dec_prt = kw_lens[1]["center_x"],kw_lens[1]["center_y"]
        X_prt,Y_prt    = conv_radec_to_xy(setting_name,ra=ra_prt,dec=dec_prt)
        if also_fits:
            str_reg+="point("+str(X_prt[0]+1)+",	"+str(Y_prt[0]+1)+") # point=x 20 color=red dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Pert. pos. from PSO } color=red\n"
            coord_y+=2        
        to_plot.append({"x":X_prt[0],"y":Y_prt[0],"rds":5,"txt":None,"col":"red","lns":"cross","lbl":"Pert. pos. from PSO"})

        
    if RM:
        kw_res = get_kwres(setting_name,updated=True)["kwargs_results"]
        ra_res,dec_res = [kw_res["ra_image_"+str(i)] for i in range(4)],[kw_res["dec_image_"+str(i)] for i in range(4)] 
        X_res,Y_res    = conv_radec_to_xy(setting_name,ra=ra_res,dec=dec_res)
        if also_fits:
            for i in range(len(X_res)):
                str_reg+="circle("+str(X_res[i]+1)+",	"+str(Y_res[i]+1)+", "+str(rd_i[i]/3.)+") # color=red dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Image pos. from MCMC } color=red\n"
            coord_y+=2
        for i in range(len(X_res)):
            rd = to_plot[i]["rds"]/3.            
            if i==0:
                lbl_im = "Image pos. from MCMC"
            else:
                lbl_im = None
            to_plot.append({"x":X_res[i],"y":Y_res[i],                 "rds":rd,"txt":None,"col":"deeppink","lns":"--",                 "lbl":lbl_im})
        
        #### add source ###
        kw_source = load_whatever(get_savefigpath(sets)+"/read_source.data")
        ra_src,dec_src = kw_source["source_ra"][0],kw_source["source_dec"][0]
        X_src,Y_src    = conv_radec_to_xy(setting_name,ra=ra_src,dec=dec_src)
        
        sig_src        = np.sqrt(np.sum(kw_source["source_ra"][1:])**2 +  np.sum(kw_source["source_dec"][1:])**2)/2.
        sig_src_xy     = sig_src/sets.pix_scale
        if also_fits:
            str_reg+="circle("+str(X_src[0]+1)+",	"+str(Y_src[0]+1)+", "+str(sig_src_xy)+") # color=chocolate dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Source pos. from MCMC } color=chocolate\n"
            coord_y+=2        
        to_plot.append({"x":X_src[0],"y":Y_src[0],"rds":sig_src_xy,"txt":None,"col":"chocolate","lns":"--","lbl":"Source pos. from MCMC"})
        
        ##### add main lens #####
        if check_if_SUB(sets):
            ra_lns,dec_lns = kw_res["center_x_lens0"],kw_res["center_y_lens0"]
        else:
            ra_lns,dec_lns = kw_res["center_x_lens_light0"],kw_res["center_y_lens_light0"]
        X_lns,Y_lns    = conv_radec_to_xy(setting_name,ra=ra_lns,dec=dec_lns)
        if also_fits:
            str_reg+="point("+str(X_lns[0]+1)+",	"+str(Y_lns[0]+1)+") # point=x 20 color=darkgray dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Main lens pos. from MCMC } color=darkgray\n"
            coord_y+=2        
        to_plot.append({"x":X_lns[0],"y":Y_lns[0],                 "rds":0,"txt":None,"col":"darkgray","lns":"cross",                 "lbl":"Main lens pos. from MCMC"})
        
        #### add vector indicating pointing angle of Lens Mass from MCMC ####
        e1,e2   = kw_res["e1_lens0"],kw_res["e2_lens0"]
        _,phi   = qphi_from_e1e2(e1,e2,ret_deg=True)
        rotang  = get_rotangle(sets) 
        phi_rt  = - phi + rotang 
        lenght  = 10 
        if also_fits:
            str_reg+="# vector("+str(X_lns[0]+1)+",	"+str(Y_lns[0]+1)+"," +str(lenght)+"," + str(phi_rt)+")  vector=1 color=cyan \n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Main lens P.A. from PSO } color=cyan \n"
            coord_y+=2
        
        ##### add perturber #####
        if check_if_SUB(sets):
            ra_prt,dec_prt = kw_res["center_x_lens_light0"],kw_res["center_y_lens_light0"]
        else:
            ra_prt,dec_prt = kw_res["center_x_lens_light1"],kw_res["center_y_lens_light1"]

        X_prt,Y_prt = conv_radec_to_xy(setting_name,ra=ra_prt,dec=dec_prt)
        if also_fits:
            str_reg+="point("+str(X_prt[0]+1)+",	"+str(Y_prt[0]+1)+") # point=x 20 color=aquamarine dash=1\n"
            str_reg+="# text("+str(int(numPix-15))+","+str(int(numPix-coord_y))+") text={Pert. pos. from MCMC } color=aquamarine\n"
            coord_y+=2        
        to_plot.append({"x":X_prt[0],"y":Y_prt[0],                 "rds":0,"txt":None,"col":"aquamarine","lns":"cross",                 "lbl":"Pert. pos. from MCMC"})
    
    if also_fits:
        str_reg+="# text(20,"+str(int(numPix-4))+") text={Values from "+strip_setting_name(setting_name)+" } font=\"12\" color=white \n"
        reg_file = get_savefigpath(sets)+"/init_img.reg"
    
        with open(reg_file, "w") as text_file:
            print(str_reg, file=text_file)

        print("Reg file in ",reg_file)
    
    image_path = sets.data_path+sets.image_name 
    fits_image = load_fits(image_path)
    fig, ax = plt.subplots(figsize=(9,9))    
    ax.imshow(fits_image, cmap = plt.cm.gist_heat)
    
    rnd = False
    for i in range(len(to_plot)):
        x,y = to_plot[i]["x"],to_plot[i]["y"]#left bottom corner, no more needed to correct for DS9 starting from 1
        rds = to_plot[i]["rds"]
        txt = to_plot[i]["txt"]
        lbl = to_plot[i]["lbl"]
        col = to_plot[i]["col"]
        lns = to_plot[i]["lns"] #linestyle
        #ax.add_patch( Rectangle((x,y),2*rds, 2*rds,color=col,fill=False,ls=lns,label=lbl) )
        if lns != "cross" and rds>.3:
            ax.add_patch( Circle((x,y),radius=rds, color=col,fill=False,ls=lns,label=lbl) )
        else:
            ax.scatter(x,y,c=col,label=lbl,marker="x")
        if txt:
            if rnd is False:
                fact=1
            else:
                fact=-1
            ax.text(x+fact*rds,y-fact*3,txt,c=col)
            rnd = not rnd
    ax.set_title("Resulting position for "+strip_setting_name(sets))
    plt.axis('off')
    plt.legend()
    plt.tight_layout()  
    if save:
        plt.savefig(get_savefigpath(sets)+"/position_image.png")#, bbox_inches='tight',pad_inches = 0)
        print("Position image in "+get_savefigpath(sets)+"/position_image.png")
    else:
        return ax


if __name__=="__main__":
    ################################
    present_program(sys.argv[0])
    ################################
    parser = argparse.ArgumentParser(description="Produce .reg file for DS9 with the initial position of images and their freedom and corresponding .png image")
    parser.add_argument("-rp", "--res_pso", dest="res_pso", default=True, action="store_false",
                        help="also plot the position resulting from the PSO")
    parser.add_argument("-rm", "--res_mcmc", dest="res_mcmc", default=True, action="store_false",
                        help="also plot the position resulting from the MCMC")

    parser.add_argument('SETTING_FILE',nargs=(1),default=[],help="setting file to consider")
    args  = parser.parse_args()
    RP    = args.res_pso
    RM    = args.res_mcmc
    setting_name =  args.SETTING_FILE[0]
    
    if RP or RM:
        from get_res import *
    create_pos_img(setting_name,RP,RM,also_fits=True,save=True)
    success(sys.argv[0])

