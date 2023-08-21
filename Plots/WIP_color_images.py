"""
The idea is to have a single program that align the images and create a color plot
problem is, for that we need resampling
    - first solution:
        - just re-define the WCS in the header such that they are compatible 
    - second solution:
        - use reproject from https://reproject.readthedocs.io/en/stable/
"""
import sys
from astropy.io import fits
from argparse import ArgumentParser

from Utils.tools import *
from Data.conversion import *
from Data.image_manipulation import get_header,shift_astrometry

@check_setting
def get_xy0(setting):
    im_path_ref = setting.data_path+setting.image_name
    x0,y0    = get_header(im_path_ref,"CRPIX1"),get_header(im_path_ref,"CRPIX2")
    return x0,y0

@check_setting
def get_radec0(setting):
    im_path_ref = setting.data_path+setting.image_name
    ra0,dec0 = get_header(im_path_ref,"CRVAL1"),get_header(im_path_ref,"CRVAL2")
    return ra0,dec0
"""
@check_setting
def correct_WCS(sett_ref,sett_list):
    sett_list = get_setting_module(sett_list,1)
    #x_ref,y_ref = sett_ref.x_image,sett_ref.y_image
    #imAx,imAy   = x_ref[0],y_ref[0]
    ra0_ref,dec0_ref = get_radec0(sett_ref)
    x0_ref,y0_ref    = get_xy0(sett_ref)
    ra_x0_ref,dec_y0_ref = conv_xy_to_radec(sett_ref,x0_ref,y0_ref)
    ra_x0_ref,dec_y0_ref = ra_x0_ref[0],dec_y0_ref[0]
    #raA,decA = conv_xy_to_radec(sett_ref,imAx,imAy)
    #raA,decA = np.array([raA,decA])+np.array([ra0,dec0])
    

    # First: find coord of A in each setting in ABS. ra,dec coord
    # by constr., sett.ra_at_xy0 and dec are the opposite of the distance in arcsec
    # of A from pix 0,0


    for sett_i in sett_list:
        if sett_i == sett_ref:
            continue
        ra0_i,dec0_i = get_radec0(sett_i)
        x0_i,y0_i    = get_xy0(sett_i)
        ra_x0,dec_y0 = conv_xy_to_radec(sett_i,x0_i,y0_i)
        ra_x0,dec_y0 = ra_x0[0],dec_y0[0] 
        dra0,ddec0     = ra0_i - ra0_ref,dec0_i - dec0_ref
        dra_x0,ddec_y0 = ra_x0 - ra_x0_ref, dec_y0 - dec_y0_ref
        Dra,Ddec       = dra0  - dra_x0   , ddec0  - ddec_y0
"""

@check_setting
def correct_WCS_simple(sett_ref,sett_list):
    # ignore real coord, set A image to ra,dec = 0,0
    sett_list = get_setting_module(sett_list,1)
    x_ref,y_ref = sett_ref.x_image,sett_ref.y_image
    imAra,imAdec   = x_ref[0],y_ref[0]
    imAx,imAy      = conv_radec_to_xy(sett_ref,imAra,imAdec)
    imAx,imAy      = imAx[0]+1,imAy[0]+1
    #raA_ref,decA_ref    = conv_xy_to_radec(sett_ref,imAx,imAy)
    #ra0_ref,dec0_ref = get_radec0(sett_ref)
    #x0_ref,y0_ref    = get_xy0(sett_ref)
    #ra_x0_ref,dec_y0_ref = conv_xy_to_radec(sett_ref,x0_ref,y0_ref)
    #ra_x0_ref,dec_y0_ref = ra_x0_ref[0],dec_y0_ref[0]
    #raA,decA = conv_xy_to_radec(sett_ref,imAx,imAy)
    #raA,decA = np.array([raA,decA])+np.array([ra0,dec0])
    
    old_name = sett_ref.data_path+sett_ref.image_name
    new_name = old_name.replace(".fits","_WCScorr_simp.fits")
    fits_ref     = copy_fits(old_name)
    fits_ref     = shift_astrometry(fits_ref,[0,0],[imAx,imAy])
    print(f"Saving {new_name}")
    fits_ref.writeto(new_name, overwrite=True)

    for sett_i in sett_list:
        if sett_i == sett_ref:
            continue
        x_i,y_i       = sett_i.x_image,sett_i.y_image
        imAx_i,imAy_i = x_i[0],y_i[0]
        imAra_i,imAdec_i   = x_i[0],y_i[0]
        imAx_i,imAy_i      = conv_radec_to_xy(sett_i,imAra_i,imAdec_i)
        imAx_i,imAy_i      = imAx_i[0]+1,imAy_i[0]+1
        old_name_i = sett_i.data_path+sett_i.image_name
        new_name_i = old_name_i.replace(".fits","_WCScorr.fits")
        fits_i     = copy_fits(old_name_i)
        fits_i     = shift_astrometry(fits_i,[0,0],[imAx_i,imAy_i])
        print(f"Saving {new_name_i}")
        fits_i.writeto(new_name_i, overwrite=True)
    print("Done")



if __name__=="__main__":
    present_program(sys.argv[0])
    parser = ArgumentParser(description="Create fits file with delensed source in colors")
    parser.add_argument('REF_SETTING_FILE',nargs=1,default="",help="Reference setting file")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    settings_list = get_setting_module(args.SETTING_FILES,sett=True)
    settings_ref  = get_setting_module(args.REF_SETTING_FILE[0],sett=True)
    correct_WCS_simple(settings_ref,sett_list=settings_list)
    success(sys.argv[0])