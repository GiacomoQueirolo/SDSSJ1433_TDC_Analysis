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
from Data.conversion import conv_xy_to_radec
from MassToLight.grid_class import Circ_Grid,Length
    
from Data.image_manipulation import fits_with_copied_hdr


# In[ ]:


def _get_EE(setting):
    setting   = get_setting_module(setting,1)
    data_path = setting.data_path
    data_path = "./data/"#"/".join(data_path.split("/")[:-1]) # main dir of all filters
    filt      = setting.filter_name.upper()
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
    flux = np.sum([lens_light_pix[i[0]][i[1]] for i in flat_grid_degraded])
    return flux

def create_PSF(sett):
    print_setting(sett)
    sett   = get_setting_module(sett,1)
    ######################
    # PSF_SS: not the best solution, but somewhat alright for now
    if float(sett.pssf)>1.:
        sett.pix_scale = sett.pix_scale/sett.pssf
        sett.transform_pix2angle = sett.transform_pix2angle / sett.pssf
    ######################
    sett.ra_at_xy0,sett.dec_at_xy0 = [0,0]

    aper,EE   = _get_EE(sett)
    radius_pix = int(Length(max(aper),sett,input_type="radec").length_pix )
    if radius_pix%2==0:
        radius_pix+=1
    new_psf    = np.zeros([radius_pix,radius_pix])

    grid = Grid_Class(center=[int(radius_pix/2),int(radius_pix/2)],radius=radius_pix+1,precision=.5,setting=sett,input_type="pixel")
    
    for i_aper,r in enumerate(aper):
        print("aperture ",r)
        flux = []
        for ditther_i in range(10):
            new_center =  np.random.normal(grid.Center.get_radec(),grid.precision_radec*.5).tolist()
            subgrid = grid.get_subgrid_with_radius_from_nc(r,new_center,input_type="radec")
            subgrid_deg_pixel =  grid._get_degraded_pixel_grid_flat(grid._flatten(subgrid))
            _flux = EE[i_aper]
            for i in subgrid_deg_pixel:
                try:
                    new_psf[i[0]][i[1]] +=_flux
                except IndexError:
                    print("subgrid larger then psf model, ignoring this point")
    new_psf/=np.sum(new_psf)
    psf_path = sett.data_path+sett.psf_name 
    new_psf_name = sett.data_path+"/test_EE_psf.fits" 
    fits_with_copied_hdr(new_psf,psf_path,data_object="",fits_res_namepath=new_psf_name,overwrite=True,verbose=True)
    
    
    plt.imshow(new_psf)
    plt.savefig(get_savefigpath(sett)+"/test_EE_psf.png")
    print("Saved "+get_savefigpath(sett)+"/test_EE_psf.png")
    plt.close()


if __name__=="__main__":
    present_program(sys.argv[0])

    parser = argparse.ArgumentParser(description="Check PSF consistency with expected value")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting files to consider")
    
    args      = parser.parse_args()
    setting_names  = args.SETTING_FILES
    if len(setting_names)>1:
        pool_obj = multiprocessing.Pool()
        _ = pool_obj.map(create_PSF,setting_names)
    else:
        create_PSF(setting_names[0])
    success(sys.argv[0])

