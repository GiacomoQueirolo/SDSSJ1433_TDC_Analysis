"""
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle(218.36022,60.10712,0.158") #text={1}
circle(218.36421,60.10829,0.158") #text={2}
circle(218.3629,60.11525,0.158") #text={3}
"""
import numpy as np
from image_manipulation import load_fitshead, get_transf_matrix
from grid_class import Point, Length, dist

def shift_grid(grid,shift,input_type="radec"):
    # shift grid center by a certain ammount
    if not isinstance(shift,Length) and type(shift,float):
        shift = Length(shift,grid.setting,input_type)
    else :
        raise ValueError("Give valid shift")
    
    
def create_ds9_circgrid(grid,file,reg_name="./ds9.reg",ref_thrsh=.05,verbose=False):
    circ = grid.circ_grid_flat
    head = load_fitshead(file)
    # center of grid in arcsec 
    center_ra,center_dec = grid.Center.ra,grid.Center.dec
    # center of grid in pix
    center_pxx,center_pxy = grid.Center.x_pix,grid.Center.y_pix
    
    ref_ra,ref_dec = head["CRVAL1"]*3600,head["CRVAL2"]*3600.
    px0_x,px0_y = head["CRPIX1"],head["CRPIX2"]
    
    # the header might have been corrected in order to have the coord referred to image A (as the lenstronomy setting is by default to do)
    # we want to check the distance in pixel from the reference point in the header and the reference point in the setting
    sett_ref_pix = np.reshape(conv_radec_to_xy(grid.setting,[0,0]),2).tolist()
    if dist([px0_x,px0_y],sett_ref_pix)<ref_thrsh:
        if verbose:
            print("The header and the setting file have the same pixel as reference point")
        shift_grid(grid,)

    else:
        transfM     = get_transf_matrix(grid.setting,in_arcsec=True) 
        ref_point   = Point(ref_ra,ref_dec)
        grid = shift_grid()

