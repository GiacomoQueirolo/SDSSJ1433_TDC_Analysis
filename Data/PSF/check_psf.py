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
import multiprocessing

from Utils.tools import *
from Data.input_data import init_kwrg_psf
from MassToLight.grid_class import Circ_Grid,Length
from Data.conversion import conv_xy_to_radec
        
from Utils.math_tools import sqrt_sum_list

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
    # we only want the unique values, so that we can count the occurence
    unique_pix = np.unique(pix_index,axis=0)
    
    resampled_psf_grid = []
    for p in unique_pix:
        # how many subpixel are there in this pixel
        hmsp = len(np.where(pix_index==p)[0])
        if hmsp>sett.pssf:
            resampled_psf_grid.append(p)
    return resampled_psf_grid
    
def compare_EE(sett):
    print_setting(sett)
    sett = get_setting_module(sett,1)
    ###########
    # PSF_SS: not the best solution, but somewhat alright for now
    if float(sett.pssf)>1.:
        sett.pix_scale = sett.pix_scale/sett.pssf
        sett.transform_pix2angle = sett.transform_pix2angle / sett.pssf
    ###########
    kwrg_psf   = init_kwrg_psf(sett,saveplots=False)
    psf_image  = kwrg_psf["kernel_point_source"]
    err_psf    = kwrg_psf["psf_error_map"]
    
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
            rnd_bound = .5
            newCenter =  np.random.normal(grid.Center.get_radec(),daper[r_i]*.5).tolist()
            while newCenter in previous_center:
                rnd_bound+=.1
                if rnd_bound>1:
                    break
                newCenter = np.random.normal(grid.Center.get_radec(),grid.edge_radec*rnd_bound).tolist()
            if rnd_bound>1:
                break
            subgrid = grid.get_subgrid(r,newCenter=newCenter,input_type="radec")
            subgrid_deg_pixel = subgrid.get_degraded_pixel_circ_grid()
            if int(sett.pssf)==1:
                Err_subgrid_deg_pixel = subgrid_deg_pixel
            else:
                Err_subgrid_deg_pixel = resample_psf_grid(sett,subgrid_deg_pixel)
            _flux = np.sum([psf_image[i[0]][i[1]] for i in subgrid_deg_pixel])
            _err  = [err_psf[i[0]][i[1]] for i in Err_subgrid_deg_pixel]
            flux.append(_flux)
            err_flux.append(sqrt_sum_list(_err))#=sqrt(SUM(sig_pixi**2))
            previous_center.append(newCenter)
        EE_psf.append(np.average(flux))
        err_EE.append(sqrt_sum_list(err_flux)/np.sqrt(len(flux)-1))
    EE_psf = np.array(EE_psf)*max(EE)/max(EE_psf)	
    err_EE = np.array(err_EE)*max(EE)/max(EE_psf)
    plt.scatter(aper,EE,label="theo. EE")
    plt.scatter(aper,EE_psf,color="darkorange",label="measured EE")
    # error is from the scatter obtained form the dither
    plt.errorbar(aper,EE_psf,yerr=err_EE,color="darkorange",fmt="x",capsize=2)
    plt.title("EE theoretic vs measured from PSF model")
    plt.xlabel("Aperture [\"]")
    plt.ylabel("Norm. Flux [pix count]")
    plt.legend()
    plt.savefig(get_savefigpath(sett)+"/EE_dither_pll.png")
    print("Saved "+get_savefigpath(sett)+"/EE_dither_pll.png")
    #plt.savefig(get_savefigpath(sett)+"/EE.png")
    #print("Saved "+get_savefigpath(sett)+"/EE.png")
    plt.close()


if __name__=="__main__":
    present_program(sys.argv[0])

    parser = argparse.ArgumentParser(description="Check PSF consistency with expected value")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting files to consider")
    
    args      = parser.parse_args()
    setting_names  = args.SETTING_FILES
    
    #for sett in setting_names:
    pool_obj = multiprocessing.Pool()
    _ = pool_obj.map(compare_EE,setting_names)
    success(sys.argv[0])
 

