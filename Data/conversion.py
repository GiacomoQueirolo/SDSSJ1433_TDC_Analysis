# if we consider the new image we convert the pix_coord of the objects
# and after that obtain the angle via pix_scale
# !! Attention: we also correct here for the difference of starting indices for array in python and ds9
import warnings
import numpy as np
from copy import copy
from astropy import units as u
from Utils.tools import get_setting_module,check_setting


#Idea: we need to change the frame of reference such that all filters has the same coordinates
# the only idea I have is to consider A as center for the F.O.R. since it's the easiest position the measure  
def conv_5 (x_im_pix,y_im_pix, \
	center_x_pix,center_y_pix,\
	x_pert_pix,y_pert_pix,\
	transform_pix2angle,ra_at_xy_0,dec_at_xy_0):
	#conversion of pixel positions in relative rotated coord
	pix_scale = ((transform_pix2angle[0][0]**2 + transform_pix2angle[0][1]**2)**.5) # "/pix
 
	#Correction for array starting at 1 in DS9 maporcavacca
	x_im_pix_conv = np.array(x_im_pix)-1
	y_im_pix_conv = np.array(y_im_pix)-1

	x_pert_pix_conv =  x_pert_pix - 1
	y_pert_pix_conv =  y_pert_pix - 1

	center_x_pix_conv = center_x_pix - 1
	center_y_pix_conv = center_y_pix - 1
	
	#conv in arcsec
	def transf(x,y,trans_mat):
		tr_x,tr_y = trans_mat.dot([x,y]).T
		return (tr_x,tr_y)
	x_image ,y_image = [],[]
	for i in range(len(x_im_pix_conv)):
		t_x,t_y = transf(x_im_pix_conv[i],y_im_pix_conv[i],transform_pix2angle)
		x_image.append(t_x+ra_at_xy_0)
		y_image.append(t_y+dec_at_xy_0)
	x_pert, y_pert =  transf(x_pert_pix_conv,y_pert_pix_conv,transform_pix2angle)
	x_pert += ra_at_xy_0
	y_pert += dec_at_xy_0
	center_x,center_y = transf(center_x_pix_conv,center_y_pix_conv,transform_pix2angle)
	center_x += ra_at_xy_0
	center_y += dec_at_xy_0
	
	ra_at_xy_0 = (x_image[0])
	dec_at_xy_0= (y_image[0])

	x_conv,y_conv = [],[]
	for i in range(len(x_image)):
		x_conv.append(x_image[i]-ra_at_xy_0)
		y_conv.append(y_image[i]-dec_at_xy_0)
	x_image,y_image=x_conv,y_conv

	x_pert -= ra_at_xy_0
	y_pert -= dec_at_xy_0

	center_x -= ra_at_xy_0
	center_y -= dec_at_xy_0
	
	ra_at_xy_0 *=-1
	dec_at_xy_0 *=-1
	return x_image, y_image, x_pert, y_pert, center_x, center_y, pix_scale, ra_at_xy_0, dec_at_xy_0	


def transf(x,y,trans_mat):
    tr_x,tr_y = trans_mat.dot([x,y]).T
    return (tr_x,tr_y)

def conv_general(ref_pos, to_convert_pos,transform_pix2angle,angle_at_pix0):
    """ 
    ref_pos = [X_ref,Y_ref] reference coord (image A)
    to_convert_pos = [[X_i,Y_i],[X_j,Y_j],...] 
    transform_pix2angle = [[C11,C12],[C21,C22]]
    angle_at_pix0 = [ra_at_xy_0,dec_at_xy_0]
    """
    #conversion of pixel positions in relative rotated coord
    #pix_scale = get_pix_scale(transform_pix2angle,only_one=True)
    pix_scale = ((transform_pix2angle[0][0]**2 + transform_pix2angle[0][1]**2)**.5)
    #Correction for array starting at 1 in DS9 maporcavacca
    to_convert_pos = np.array([np.array([k-1. for k in i]) for i in to_convert_pos])
    ref_pos        = np.array([ref_pos[0]-1,ref_pos[1]-1])
    new_to_conv    = np.zeros_like(to_convert_pos)
    
    #conv in arcsec
    for i in range(len(to_convert_pos)):
        t_x,t_y        = transf(*to_convert_pos[i],transform_pix2angle)
        new_to_conv[i] = [t_x+angle_at_pix0[0],t_y+angle_at_pix0[1]]
    
    ref_conv = [transf(*ref_pos,transform_pix2angle)[i]+angle_at_pix0[i] for i in range(len(ref_pos))]
    angle_at_pix0 = copy(np.array(ref_conv))
    for i in range(len(new_to_conv)):
        new_to_conv[i] = new_to_conv[i] - angle_at_pix0
    return new_to_conv,pix_scale,-1*angle_at_pix0

def conv_general_easy(ref_pos, to_convert_pos,transform_pix2angle,angle_at_pix0,Aap0_in_deg=True):
    """ 
    ref_pos = [X_ref,Y_ref] reference coord (image A)
    to_convert_pos = [[X_i,Y_i],[X_j,Y_j],...] 
    transform_pix2angle = [[C11,C12],[C21,C22]]
    angle_at_pix0 = [ra_at_xy_0,dec_at_xy_0] 
    Aap0_in_deg   = if angle_at_pix0 is in degrees or arcsecs
    """
    #conversion of pixel positions in relative rotated coord
    #pix_scale = get_pix_scale(transform_pix2angle,only_one=True)
    pix_scale = ((transform_pix2angle[0][0]**2 + transform_pix2angle[0][1]**2)**.5)
    #Correction for array starting at 1 in DS9 maporcavacca
    to_convert_pos = np.array([np.array([k-1. for k in i]) for i in to_convert_pos])
    ref_pos        = np.array([ref_pos[0]-1,ref_pos[1]-1])
    new_to_conv    = np.zeros_like(to_convert_pos)
    scale          = 1
    if Aap0_in_deg:
        scale=3600
    #conv in arcsec
    for i in range(len(to_convert_pos)):
        t_x,t_y        = transf(*to_convert_pos[i],transform_pix2angle)
        new_to_conv[i] = np.array([t_x,t_y])+scale*np.array(angle_at_pix0)
    return new_to_conv,pix_scale
"""
@check_setting
def conv_total(setting,to_conv_pos):
    transform_pix2angle = get_transf_matrix(setting.data_path+setting.image_name,in_arcsec=True)
    angle_at_pix0 = [0,0] # should be always 0,0
    return conv_general(to_conv_pos,transform_pix2angle,angle_at_pix0)
"""

def find_center(x_im,y_im): # old
    # give pix coord of images
    # return pix coord of the "center" of the image
    X = .5*(max(x_im)+min(x_im))
    Y = .5*(max(y_im)+min(y_im))
    return X,Y

def find_barycenter(x_im,y_im):
    # give pix coord of images
    # return pix coord of the barycenter of the image
    X = sum(x_im)/len(x_im)
    Y = sum(y_im)/len(y_im)
    return X,Y

    
###### Invert conversion #################

def inv_transf(ra,dec,trans_mat):
    inv_trans_mat = np.linalg.inv(trans_mat)
    x,y = inv_trans_mat.dot([ra,dec]).T
    return (x,y)

def invert_conv(ra_im,dec_im,A_ra,A_dec,transform_pix2angle):    
    ra_im = np.array(ra_im)+A_ra
    dec_im = np.array(dec_im)+A_dec
    
    x_im =[]
    y_im =[]
    try:
        for i in range(len(ra_im)):
            tr_x,tr_y = inv_transf(ra_im[i],dec_im[i],transform_pix2angle)
            x_im.append(tr_x+1)
            y_im.append(tr_y+1)
    except:
        tr_x,tr_y = inv_transf(ra_im,dec_im,transform_pix2angle)
        x_im.append(tr_x+1)
        y_im.append(tr_y+1)
    # this are in DS9 frame of reference (+1 wrt python)
    return x_im,y_im

@check_setting
def conv_radec_to_xy(setting,ra,dec,ref_radec=None):
    try:
        ra,dec = np.array(ra),np.array(dec) 
    except:
        ra,dec = np.array([ra]),np.array([dec]) 
    # by construction ra/dec_at_xy_0 are defined by pos of A image, which is now in the 0,0 coord
    if ref_radec is None:
        neg_ra_ref  = -setting.ra_at_xy_0
        neg_dec_ref = -setting.dec_at_xy_0
    else:
        neg_ra_ref  = -ref_radec[0]
        neg_dec_ref = -ref_radec[1]
    X,Y = invert_conv(ra,dec,neg_ra_ref,neg_dec_ref,setting.transform_pix2angle)
    # convert to python coord
    X,Y = np.array(X)-1,np.array(Y)-1
    return X.tolist(),Y.tolist()
"""
def conv_xy_to_radec(setting,x,y):
    sets      = get_setting_module(setting,sett=True)
    Ara,Adec  = sets.ps_params[0][0]["ra_image"][0],sets.ps_params[0][0]["dec_image"][0]
    Axy       = np.reshape(conv_radec_to_xy(sets,Ara,Adec),2).tolist()
    radec,_,_ = conv_general(Axy, [[x,y]] ,sets.transform_pix2angle,[sets.ra_at_xy_0,sets.dec_at_xy_0])
    return radec
"""
@check_setting
def conv_xy_to_radec(setting,x,y,xy_ref=None):
    # x,y are given in DS9 coordinates
    if xy_ref is None:
        # this should be 0 by construction:
        #Ara,Adec  = sets.ps_params[0][0]["ra_image"][0],sets.ps_params[0][0]["dec_image"][0]
        #xy_ref    = np.reshape(conv_radec_to_xy(sets,Ara,Adec),2)
        # the /2 is bc it is summed with himself if not (same res. with 2* the input and *-1 the output)
        xy_ref = np.reshape(np.array(conv_radec_to_xy(setting,ra=-setting.ra_at_xy_0,dec=-setting.dec_at_xy_0))/2.,2)
    radec,_,_ = conv_general(xy_ref, [[x,y]] ,setting.transform_pix2angle,[setting.ra_at_xy_0,setting.dec_at_xy_0])
    ra,dec    = [radec_i[0] for radec_i in radec],[radec_i[1] for radec_i in radec]
    return ra,dec
    
def conv_xy_from_setting1_to_setting2(x1,y1,setting1,setting2):
    """
    x1,y1    = pixel coordinatex with respect to setting1
    setting1 = name of setting file for coord x1,y1
    setting2 = name of setting file over which we want to convert
    """
    ra,dec = conv_total(setting1,[x1,y1])
    x2,y2  = conv_radec_to_xy(setting2,ra,dec)
    return x2,y2

### from ellipticities_conversion.py -> heaviliy reworked from param_util.py from lenstronomy ###

def e1e2_from_qphi(q,phi,deg=True):
    """
    transforms orientation angle and axis ratio into complex ellipticity moduli e1, e2

    :param phi: angle of orientation (in radian)
    :param q: axis ratio minor axis / major axis
    :return: eccentricities e1 and e2 in complex ellipticity moduli
    """
    if deg:
        warnings.warn("NOTE: assumed phi given in degree!")
        phi = np.array(phi)*u.deg.to("rad")

    c  = (1. - q) / (1. + q) 
    e1 = np.cos(2 * phi)*c
    e2 = np.sin(2 * phi)*c
    return e1, e2


def qphi_from_e1e2(e1, e2,ret_deg=False):
    """
    transforms complex ellipticity moduli in orientation angle and axis ratio

    :param e1: eccentricity in x-direction
    :param e2: eccentricity in xy-direction
    :return: angle in radian, axis ratio (minor/major)
    """
    phi = np.arctan2(e2, e1)/2
    c = np.sqrt(e1**2+e2**2)
    #c = np.minimum(c, 0.9999)
    q = (1-c)/(1+c)
    #q = (1-np.sqrt(e1**2+ e2**2))/(1+np.sqrt(e1**2+ e2**2))
    if ret_deg:
        warnings.warn("NOTE: phi is returned in degrees!")
        return q,phi*u.rad.to("deg")
    else:
        warnings.warn("NOTE: phi is returned in radians!")
        return q,phi    
#########################
# MOD_VISPOS
########################
def get_pixscale(transform_pix2angle):
    pix_scale = ((transform_pix2angle[0][0]**2 + transform_pix2angle[0][1]**2)**.5) # "/pix
    return pix_scale

def conv_std(vec_pix,Apixcoord,transform_pix2angle):
    # vec_pix: vector of pixel coordinates, (N,2) in DS9 coordinates
    # Apixcoord: vector of pixel coordinates of image A for this image, (2) in DS9 coordinates
    # standard pos conversion 
    #Correction for array starting at 1 in DS9 maporcavacca
    vec_pix   = np.array(vec_pix)-1
    Apixcoord = np.array(Apixcoord)-1
    pix_wrt_A = vec_pix - Apixcoord
    
    vec  = pix_wrt_A.dot(np.array(transform_pix2angle).T) 
    # to obtain radec_at_xy0 imput 0,0
    return vec

def conv_VISPOS(vec_pix,Apixcoord,transform_pix2angle,R_VISPOS,t_VISPOS):
    # vec_pix: vector of pixel coordinates, (N,2) in DS9 coordinates
    # Apixcoord: vector of pixel coordinates of image A for this image, (2) in DS9 coordinates
    # transform_pix2angle: the transformation matrix for this image
    # R_VISPOS,t_VISPOS: the rotation and the translation vector to match this image to the VIS reference frame 
    
    #Correction for array starting at 1 in DS9 maporcavacca
    vec_pix = np.array(vec_pix)-1
    Apixcoord = np.array(Apixcoord)-1
    pix_wrt_A = vec_pix - Apixcoord
    
    # combined transformation matrix (already transposed)
    vec_transl = np.array(t_VISPOS).dot(np.array(R_VISPOS).T)
    MT   = np.array(transform_pix2angle).T.dot(np.array(R_VISPOS).T)
    vec  = pix_wrt_A.dot(MT) + vec_transl

    # to obtain radec_at_xy0 input 1,1
    return vec

# remember, the pix coord all start from 1,1 everytime

def conv_sett(sett,vec_pix):
    vec_radec_at_0 = sett.ra_at_xy_0,sett.dec_at_xy_0
    transform_pix2angle = sett.transform_pix2angle
    Apixcoord  = invconv_std([0,0],vec_radec_at_0,transform_pix2angle)
    return conv_std(vec_pix,Apixcoord,transform_pix2angle)

def invconv_std(vec_radec,radec_at_xy0,transform_pix2angle):
    # vec_radec: vector in radec coordinates, (N,2), wrt A
    # radec_at_xy0: vector in radec coordinates of image pix 0,0 for this image, (2) in DS9 coordinates
    # inverted  pos conversion 
    inv_trans_mtx = np.linalg.inv(transform_pix2angle)
    Apixcoord = -inv_trans_mtx.dot(radec_at_xy0).T 
    pix_wrt_A = inv_trans_mtx.dot(vec_radec).T
    vec_pix   = pix_wrt_A + Apixcoord 
    #Correction for array starting at 1 in DS9 maporcavacca
    vec_pix   = np.array(vec_pix)+1
    return vec_pix


def invconv_sett(sett,vec_radec):
    radec_at_xy0 = sett.ra_at_xy_0,sett.dec_at_xy_0
    transform_pix2angle = sett.transform_pix2angle
    
    return invconv_std(vec_radec,radec_at_xy0,transform_pix2angle)


def get_transf_matrix(image_path,in_arcsec=True):
    with fits.open(image_path) as hdulist:
        hdr = hdulist[0].header
    CD1_1,CD1_2,CD2_1,CD2_2 = hdr["CD1_1"],hdr["CD1_2"],hdr["CD2_1"],hdr["CD2_2"]
    transform_pix2angle = np.array([[CD1_1, CD1_2], [CD2_1, CD2_2]])
    if in_arcsec:
        transform_pix2angle*=3600.
    return transform_pix2angle
    