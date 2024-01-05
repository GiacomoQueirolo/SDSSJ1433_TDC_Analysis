#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Check that psf behaves almost as the expected one from 
#https://iopscience.iop.org/article/10.1088/0067-0049/214/2/24/pdf
# and the EE
# https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/ir-encircled-energy
# https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/uvis-encircled-energy
import sys,csv
import argparse
import numpy as np
import matplotlib.pyplot as plt


from Utils.tools import *
from Utils.get_res import load_whatever

from MassToLight.grid_class import Circ_Grid,Length
from Data.conversion import conv_xy_to_radec,conv_radec_to_xy
        
from Data.image_manipulation import load_fits

def _get_EE(setting):
    setting   = get_setting_module(setting,1)
    data_path = setting.data_path
    data_path = "./data/"#"/".join(data_path.split("/")[:-1]) # main dir of all filters
    filt      = setting.filter_name.upper()
    ####
    # correct f160w_7030
    if filt=="F160W_7030":
        filt = "F160W"
    ####
    if "F1"==filt[:2]:
        # IR
        EE_file = "ir_ee_corrections.csv"
    else:
        # UV
        EE_file = "wfc3uvis2_aper_007_syn.csv"
    EE_file = data_path+"/"+EE_file
    with open(EE_file,newline="") as f:
        reader = csv.reader(f)
        for i,EErow in enumerate(reader):
            if i==0:
                titles = EErow
                # APER is in arcsec
            if filt in EErow[0]:
                EErow_i = EErow[2:]
                break
    aper = np.array([float(aper.split("#")[1]) for aper in titles[2:]])
    EE   = np.array([float(ee) for ee in EErow_i])
    return aper,EE

def flux_in_grid(lens_light_pix,grid): 
    # now the problem is, the grid is not discretised in pixels
    # we now have a grid with the coord of the pixels
    flux = np.sum([lens_light_pix[i[0]][i[1]] for i in grid.get_degraded_pixel_grid()])
    return flux
    
def resample_psf_grid(sett,subgrid_deg_pixel):
    # we want the resampled version of the grid given the original pixel resolution
    # which we have in the error frame of the PSF
    
    # exact pixel index in the orig resolution
    index = np.array(subgrid_deg_pixel)/sett.pssf
    # integer version of the pixel index
    pix_index = index.astype(int)
    #return pix_index
    
    # we only want the unique values, so that we can count the occurence
    unique_pix = np.unique(pix_index,axis=0)
    
    resampled_psf_grid = []
    for p in unique_pix:
        # how many subpixel are there in this pixel
        hmsp = len(np.where(pix_index==p)[0])
        if hmsp>sett.pssf:
            resampled_psf_grid.append(p)
    return np.array(resampled_psf_grid)-1
    
def acceptable_newcenter(newCenter_xy,naxis,ri):
    # check that the new center is not too close to the edges
    # assuming naxis1==naxis2
    xc,yc = newCenter_xy
    if xc-ri<=0:
        return False
    if xc+ri>=naxis:
        return False
    if yc-ri<=0:
        return False
    if yc+ri>=naxis:
        return False
    return True
 
aper_name   = "psf_aper.json"
EE_name     = "/psf_theoEE.json"
EE_psf_name = "/psf_measEE.json"

def compare_EE(sett,savedata=True):
    print_setting(sett)
    sett = get_setting_module(sett,1)
    ###########
    # PSF_SS: not the best solution, but somewhat alright for now
    if float(sett.pssf)>1.:
        sett.pix_scale = sett.pix_scale/sett.pssf
        sett.transform_pix2angle = sett.transform_pix2angle / sett.pssf
    ###########
    #kwrg_psf   = init_kwrg_psf(sett,saveplots=False)
    #psf_image  = kwrg_psf["kernel_point_source"]
    psf_file  = sett.data_path+sett.psf_name 
    psf_image = np.array(load_fits(psf_file))
    
    
    aper,EE    = _get_EE(sett)
    naxis1,naxis2 = len(psf_image),len(psf_image[0])
    if naxis1!=naxis2:
        raise RuntimeError("psf is not a square?")
    naxis = naxis1
    brightest     = np.array(np.where(psf_image==np.max(psf_image))).reshape(2).tolist()
    print("center image  |   brightest pixel")
    print(naxis/2.,"  |   ",brightest)
    
    center       = [naxis/2.,naxis/2.]
    center_radec =  np.reshape(conv_xy_to_radec(sett,*center),2).tolist()
    sett.ra_at_xy0,sett.dec_at_xy0 = center_radec
    ############
    # PSF_small: we cut out the EE which are not in the PSF
    radius_pix = Length(max(aper),sett,input_type="radec").length_pix 
    if radius_pix>naxis:
        aper_pix = np.array([Length(a,sett,input_type="radec").length_pix for a in aper])
        index_EE = np.where((aper_pix<=min(center)) & (aper_pix<=(naxis-max(center))))
        aper     = np.array(aper)[index_EE]
        EE       = np.array(EE)[index_EE]
        radius_pix = Length(max(aper),sett,input_type="radec").length_pix 
    ############
    grid = Circ_Grid(center=center,radius=radius_pix+1,edge=.5,setting=sett,input_type="pixel")
    EE_psf = []
    err_EE = []
    daper = [aper[i+1]-aper[i] for i in range(len(aper[:-1]))]
    daper.append(daper[-1])
    for r_i,r in enumerate(aper):
        print("aperture ",r)
        flux = []
        err_flux = []
        previous_center=[]
        for dither_i in range(10):
            if dither_i!=0:
                rnd_bound = .5
                newCenter =  np.random.normal(grid.Center.get_radec(),np.mean(daper)*.25).tolist()
                while ( not acceptable_newcenter(conv_radec_to_xy(sett,*newCenter),naxis,r) ) or ( newCenter in previous_center):
                    rnd_bound+=.1
                    if rnd_bound>1:
                        break
                    newCenter = np.random.normal(grid.Center.get_radec(),grid.edge_radec*rnd_bound).tolist()
                if rnd_bound>1:
                    break
            else:
                newCenter = grid.Center.get_radec()
            subgrid = grid.get_subgrid(r,newCenter=newCenter,input_type="radec")
            subgrid_deg_pixel = subgrid.get_degraded_pixel_circ_grid()
            """if int(sett.pssf)==1:
                Err_subgrid_deg_pixel = subgrid_deg_pixel
            else:
                Err_subgrid_deg_pixel = resample_psf_grid(sett,subgrid_deg_pixel)"""
            try:
                _flux = np.sum([psf_image[i[0]][i[1]] for i in subgrid_deg_pixel])
                flux.append(_flux)
                previous_center.append(newCenter)
            except IndexError:
                print(f"Dither {dither_i} failed due to subgrid_pixel being out of bounds")
                continue
            #_err  = [err_psf[i[0]][i[1]]/(sett.pssf**2) for i in Err_subgrid_deg_pixel]
            #err_flux.append(sqrt_sum_list(_err))#=sqrt(SUM(sig_pixi**2))
        EE_psf.append(np.average(flux))
        #err_EE.append(sqrt_sum_list(err_flux)/np.sqrt(len(flux)-1))
    EE_psf = np.array(EE_psf)*max(EE)/max(EE_psf)
    #err_EE = np.array(err_EE)*max(EE)/max(EE_psf)
    if savedata:
        svpth=get_savemcmcpath(sett)
        print(svpth)
        print(svpth+"/"+aper_name)
        save_json(aper,svpth+"/"+aper_name)
        print(svpth+"/"+EE_name)
        save_json(EE,svpth+"/"+EE_name)
        print(svpth+"/"+EE_psf_name)
        save_json(EE_psf,svpth+"/"+EE_psf_name)
    return [aper,EE,EE_psf]
    
def plot_EE(sett,aper,EE,EE_psf):
    fnt = 16
    plt.rcParams['xtick.labelsize'] = fnt
    plt.rcParams['ytick.labelsize'] = fnt 
    plt.rcParams['font.size'] = fnt
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('font', size=fnt)          # controls default text sizes
    plt.rc('axes', titlesize=fnt)     # fontsize of the axes title
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('legend', fontsize=fnt)    # legend fontsize

    plt.figure(figsize=(6,5))
    plt.plot(aper,EE,label="theoretical EE")
    plt.scatter(aper,EE_psf,color="darkorange",label="measured EE")
    plt.title("Encircled Energy of PSF model")
    plt.xlabel("Aperture [\"]")
    plt.ylabel("Normalised Flux [pix count]")
    plt.legend(loc='lower right')
    #plt.tight_layout()
    plt.savefig(get_savefigpath(sett)+"/EE_dither_pll.pdf")
    print("Saved "+get_savefigpath(sett)+"/EE_dither_pll.pdf")
    plt.close()


if __name__=="__main__":
    present_program(sys.argv[0])

    parser = argparse.ArgumentParser(description="Check PSF consistency with expected value")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting files to consider")
    
    args      = parser.parse_args()
    setting_names  = args.SETTING_FILES
    
    for sett in setting_names:
        sett = get_setting_module(sett,1)
        try:
            svpth  = get_savemcmcpath(sett)
            aper   = load_whatever(svpth+"/"+aper_name)
            EE     = load_whatever(svpth+"/"+EE_name)
            EE_psf = load_whatever(svpth+"/"+EE_psf_name)
            print("loaded previous data")
        except:
            #pool_obj = multiprocessing.Pool()
            list_EE = compare_EE(sett)#pool_obj.map(compare_EE,sett)
            aper,EE,EE_psf = list_EE
        plot_EE(sett,aper,EE,EE_psf)    
    success(sys.argv[0])
 

