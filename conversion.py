# if we consider the new image we convert the pix_coord of the objects
# and after that obtain the angle via pix_scale
# !! Attention: we also correct here for the difference of starting indices for array in python and ds9
import numpy as np
from tools import get_setting_module
from image_manipulation import get_transf_matrix
from copy import copy
from astropy import units as u
import warnings
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
    def transf(xy,trans_mat):
        tr_x,tr_y = np.array(trans_mat).dot(xy).T
        return np.array([tr_x,tr_y])
    for i in range(len(to_convert_pos)):
        t_x,t_y        = transf(to_convert_pos[i],transform_pix2angle)
        new_to_conv[i] = [t_x+angle_at_pix0[0],t_y+angle_at_pix0[1]]
    
    ref_conv = [transf(ref_pos,transform_pix2angle)[i]+angle_at_pix0[i] for i in range(len(ref_pos))]
    angle_at_pix0 = copy(np.array(ref_conv))
    for i in range(len(new_to_conv)):
        new_to_conv[i] = new_to_conv[i] - angle_at_pix0
    return new_to_conv,pix_scale,-1*angle_at_pix0



def conv_total(setting,to_conv_pos):
    sets = get_setting_module(setting).setting()
    transform_pix2angle = get_transf_matrix(sets.data_path+sets.image_path,in_arcsec=True)
    angle_at_pix0 = [0,0] # should be always 0,0
    return conv_general(to_conv_pos,transform_pix2angle,angle_at_pix0)
    
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

def transf(x,y,trans_mat):
    tr_x,tr_y = trans_mat.dot([x,y]).T
    return (tr_x,tr_y)

def inv_transf(ra,dec,trans_mat):
    inv_trans_mat = np.linalg.inv(trans_mat)
    x,y = inv_trans_mat.dot([ra,dec]).T
    return (x,y)

def invert_conv_5 (ra_im, dec_im,A_x,A_y,transform_pix2angle): 

    ra_A, dec_A = transf(A_x-1,A_y-1,transform_pix2angle) 
    
    ra_im = np.array(ra_im)+ra_A
    dec_im = np.array(dec_im)+dec_A
    
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
    
    return x_im,y_im



def conv_radec_to_xy(setting_name,ra,dec):
    sets          = get_setting_module(setting_name).setting()
    Ara,Adec      = sets.ps_params[0][0]["ra_image"][0],sets.ps_params[0][0]["dec_image"][0]
    neg_Ax,neg_Ay = invert_conv_5([sets.ra_at_xy_0],[sets.dec_at_xy_0],Ara-1,Adec-1,sets.transform_pix2angle)
    try:
        ra,dec = np.array(ra),np.array(dec) 
    except:
        ra,dec = np.array([ra]),np.array([dec]) 
    
    X,Y = invert_conv_5(ra,dec,-neg_Ax[0],-neg_Ay[0],sets.transform_pix2angle)
    
    return X,Y

def conv_xy_from_setting1_to_setting2(x1,y1,setting1,setting2):
    """
    x1,y1    = pixel coordinatex with respect to setting1
    setting1 = name of setting file for coord x1,y1
    setting2 = name of setting file over which we want to convert
    """
    ra,dec = conv_total(setting1,[x1,y1])
    x2,y2 = conv_radec_to_xy(setting2,ra,dec)
    return x2,y2



### from ellipticities_conversion.py -> heaviliy reworked from param_util.py from lenstronomy ###

def e12_from_qphi(phi, q,deg=True):
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
