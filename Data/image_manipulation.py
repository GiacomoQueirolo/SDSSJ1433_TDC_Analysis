#!/usr/bin/env python
# coding: utf-8

#2]:


import csv
import warnings
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from Utils.tools import *
    
def mask_in(coord_x,coord_y,rad,x_grid,y_grid=None):
    if y_grid ==None:
        y_grid = x_grid
    mask = np.ones(shape=(x_grid,y_grid))
    for i in range(len(mask)):
        for j in range(len(mask[0])):
            R = np.sqrt((j-coord_x)**2 +(i-coord_y)**2)
            if R<=rad:
                mask[i][j]=0
    return mask

def mask_out(coord_x,coord_y,rad,x_grid,y_grid=None):
    if y_grid ==None:
        y_grid = x_grid
    mask = np.ones(shape=(x_grid,y_grid))
    for i in range(len(mask)):
        for j in range(len(mask[0])):
            R = np.sqrt((j-coord_x)**2 +(i-coord_y)**2)
            if R>rad:
                mask[i][j]=0
    return mask

def create_mask(image,setting):
    numPix  = image.shape[0]
    setting = get_setting_module(setting).setting()
    if not hasattr(setting,"R_tot_mask"):
        R_tot_mask = numPix*1./2.
    else:
        R_tot_mask = setting.R_tot_mask
    # we mask all the objects
    mask = np.ones(shape=(numPix,numPix))
    for i in range(len(setting.x_mask)):
        mask*=mask_in(setting.x_mask[i] ,setting.y_mask[i],setting.r_mask[i],numPix) 
    
    # we mask everything further out than a certain radius
    if not hasattr(setting,"new_center_mask_out"):
        x_mask_out,y_mask_out = numPix/2,numPix/2
    else:
        x_mask_out,y_mask_out = setting.new_center_mask_out
        
    mask*=mask_out(x_mask_out,y_mask_out,rad=R_tot_mask,x_grid=numPix)
    
    # We mask eventual residuals in center (due to lens light modelling imperfections)
    for i in range(len(setting.x_mask_cen)):
        mask*=mask_in(coord_x=setting.x_mask_cen[i],\
                      coord_y=setting.y_mask_cen[i],\
                      rad=setting.rad_mask_cen[i],\
                      x_grid=numPix)
    return mask



def load_fits(image_path,HDU=0):
    #load the image and read it as numpy array
    with fits.open(image_path) as hdulist:
        image   = hdulist[HDU].data
    return image

def load_fitshead(image_path,HDU=0):
    #load the image header
    with fits.open(image_path) as hdulist:
        head   = hdulist[HDU].header
    return head
    
def get_lens_light_model(setting):
    setting  = get_setting_module(setting).setting()
    lens_light_model_path = setting.data_path+"/"+setting.lens_light_model_name
    image_model = load_fits(lens_light_model_path,HDU=0) 
    return image_model
    
    
def match_image_sub(image,image_model,setting): 
    # WIP
    setting = get_setting_module(setting).setting()
    image_path = setting.data_path+"/"+setting.image_name
    with fits.open(image_path,ignore_missing_end=True) as hdulist:
        hdr = hdulist[0].header
    lens_light_model_path = setting.data_path+"/"+setting.lens_light_model_name
    with fits.open(lens_light_model_path,ignore_missing_end=True) as hdulist:
        model_hdr  = hdulist[0].header
    if setting.already_sub:
        history = hdr["HISTORY"]

        tocrop = []
        sub_ = False # we crop only the crops done after the subtraction
        for h in history:
            if "subtracted" in h:
                sub_ = True
            if "extracted" in h and sub_:
                tocrop.append(h.split(" ")[-1])
        try:
            history_model = model_hdr["HISTORY"]
            already_cropped = [] 
            for h in history_model:
                if "extracted" in h:
                    already_cropped.append(h.split(" ")[-1])
            if already_cropped!=[]:
                # pragma no cover
                raise             
        except KeyError:
            # no history of the model, no cropped region
            pass
        
        for tc in tocrop:       
            reg = tc.replace(",",":").split(":") # [int(r[2])-1:int(r[3]),int(r[0])-1:int(r[1])]
            reg = int(reg[0]),int(reg[2]),int(reg[1]),int(reg[3])
            image_model = image_model[reg[2]-1:reg[3],reg[0]-1:reg[1]]# from extract_fits.py of Matthias Kluge
        if any(np.shape(image_model)==0):
            raise RuntimeError("Image model not compatible, something went wrong")
    else:

        try:
            dx = model_hdr["CRPIX1"] - hdr["CRPIX1"] + 1
            dy = model_hdr["CRPIX2"] - hdr["CRPIX2"] + 1
            # pragma no cover
            raise
        except:
            print("Warning: Header key CRPIX1 or CRPIX2 not found.")
            # pragma no cover
            raise
    if np.shape(image_model)==np.shape(image):
        return image_model
    else:
        raise RuntimeError("Image model not compatible, something went wrong:",np.shape(image_model),"!=",np.shape(image)) 
               
def subtract_light(image,setting):
    setting  = get_setting_module(setting).setting()
    # Input the image and the setting file
    # Output the image with subtracted lens light
    if setting.sub==True:
        #if setting.lens_light_model_name is None:
            #print("Warning: no lens light model to subtract. I assume the input image was already subtracted.\n")
        if setting.already_sub is True:
            print("Lens light model already subtracted.\n")
            image_sub = image
        else:
            image_model = get_lens_light_model(setting)
            if not np.shape(image_model)==np.shape(image): 
                image_model = match_image_sub(image,image_model,setting)
            #We substract from the original image
            image_sub = image - image_model
    else:
        return image
    return image_sub

def get_transf_matrix(image_path_or_setting,in_arcsec=True):
    try:
        setting = get_setting_module(image_path_or_setting).setting()
        image_path = setting.data_path+"/"+setting.image_name
    except RuntimeError:    
        image_path = image_path_or_setting

    with fits.open(image_path) as hdulist:
        hdr = hdulist[0].header
    CD1_1,CD1_2,CD2_1,CD2_2 = hdr["CD1_1"],hdr["CD1_2"],hdr["CD2_1"],hdr["CD2_2"]
    transform_pix2angle = np.array([[CD1_1, CD1_2], [CD2_1, CD2_2]])
    if in_arcsec:
        transform_pix2angle*=3600.
    return transform_pix2angle


#4]:


def plot_image(image,setting,savefig_path_name=None,err_image=False):
    cmap_string = 'gist_heat'
    cmap = plt.get_cmap(cmap_string)
    cmap.set_bad(color='k', alpha=1.)
    cmap.set_under('k')
    setting  = get_setting_module(setting).setting()
    if not err_image:
        if hasattr(setting, 'v_min'):
            v_min = setting.v_min
            v_max = setting.v_max
        else:
            v_min = -4
            v_max = 1
    else:
        if hasattr(setting, 'e_v_min'):
            v_min = setting.e_v_min
            v_max = setting.e_v_max
        else:
            v_min = -1
            v_max = 0
    f, ax = plt.subplots(1, 1, figsize=(6, 6), sharex=False, sharey=False)
    ax.matshow(np.log10(image), origin='lower',  vmin=v_min, vmax=v_max, cmap=cmap, extent=[0, 1, 0, 1])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.autoscale(False)
    if savefig_path_name is None:
        return ax
    else:
        plt.savefig(savefig_path_name)

def plot_projection(data,savefigpath=None,ax0=None,udm="<e-/sec>",title="PSF profile",filename="PSF_projection",plot_midax=True):
    # by default this is implemented for the PSF projection
    # data = NxM matrix 
    # ax0  = axis 0 along which project (by default, center)
    # udm  = label of y axis
    # title = title of the plot
    # filename = filename to save
    # savefigpath = path where to save
    cnt_data_x = int(len(data[1])/2.)
    cnt_data_y = ax0 if ax0 else int(len(data[0])/2.)
    dx  = int(cnt_data_x/2)
    x0  = int(cnt_data_x-dx)
    x1  = int(cnt_data_x+dx)

    proj_data  = data[cnt_data_y][x0:x1]
    x_projdata = np.arange(x0,x1,dtype=int)
    fig,ax = plt.subplots(figsize=(6,6),dpi=80)
    ax.plot(x_projdata,proj_data,c="k")
    #,marker="x",color="k")
    ax.set_xlabel("Pixels[]")
    ax.set_ylabel(udm)
    ax.set_title(title)
    if plot_midax:
        ax.axvline(len(data)/2.,label="Image Center",ls="--",c="r")
        for max_i in x_projdata[np.array(proj_data)==np.max(proj_data)]:
            ax.axvline( float(max_i),label="Maximum pixel",ls="-.",c="b")
        ax.legend()
    if not savefigpath is None:
        plt.savefig(savefigpath+"/"+filename.replace(".pdf","")+".pdf")
        return 0
    else:
        return ax

def crop_egdes(image,setting):
    setting  = get_setting_module(setting).setting()
    # Crop image
    if not hasattr(setting, 'edge'):
        raise RuntimeError("Setting file has no edge")
    else:
        edge = getattr(setting, "edge")
    if hasattr(setting, 'dx'):
        dx = setting.dx
        dy = setting.dy
    else:
        dx = int(setting.center_x_pix)
        dy = int(setting.center_y_pix)

    x_min,x_max = dx - edge, dx + edge
    y_min,y_max = dy - edge, dy + edge
    
    image_cut = []
    for i in range(len(image)):
        im_line =[]
        if i>y_min and i<y_max:
            for j in range(len(image[i])):
                if j>x_min and j<x_max:
                    im_line.append(image[i][j])
            image_cut.append(im_line)
    image_cut = np.array(image_cut)
    return image_cut

def psf_correction(psf,setting):
    corr_psf = psf
    setting  = get_setting_module(setting).setting()
    if setting.pssf>1:
        dl_mean0 = np.mean(psf[0])
        dl_mean1 = np.mean(psf[:,0])
        corr_psf = np.ones((len(psf)+2,len(psf[0])+2))
        corr_psf[0]*=dl_mean0
        corr_psf[-1]*=dl_mean0
        corr_psf[:,0]*=dl_mean1
        corr_psf[:,-1]*=dl_mean1
        corr_psf[0][0]=dl_mean0
        corr_psf[0][-1]=dl_mean0
        corr_psf[-1][-1]=dl_mean1
        corr_psf[-1][0]=dl_mean1

        for i in range(1,len(corr_psf)-1):
            for j in range(1,len(corr_psf[0])-1):
                corr_psf[i][j]=psf[i-1][j-1]
    return corr_psf



#########################
def _acceptable(matrix_shape,coord):
    x,y = coord
    if x>=0 and y>=0:
        if x<=matrix_shape[0]-1 and y<=matrix_shape[1]-1:
            return True
    return False

class Px_Contour():
    def __init__(self,row,col,fam=None,type="mask"):
        self.row = row
        self.col = col
        self.fam = fam
        self.type = type
        self.bounds = []
    def coord(self):
        return [self.row,self.col]
    def get_bordering_coord(self):
        bord_coords = [[self.row-1,self.col],
                       [self.row,self.col-1],
                       [self.row,self.col+1],
                       [self.row+1,self.col]]
        return bord_coords
    def get_bordering_pixels(self):
        bord_coords = self.get_bordering_coord()
        boord_pxs = [Px_Contour(col=bc[0],row=bc[1],fam=self.fam) for bc in bord_coords]
        return boord_pxs
    def __eq__(self, __o: object):
        return __o.row==self.row and __o.col==self.col

def get_neighboring_coords(matrix, coords):
    # find the neighboring pixels to the given ones 
    # and to which "family" of pixel it belongs to
    # divided into one for the pixels and one for the neighbors
    #neighbors = [] #the contours pixels
    #mpix_fam  = [] #family of the single pixel 
    #neigh_fam = [] #family of the contour pixel
    index_fam   = np.arange(len(matrix)*len(matrix[1])).reshape(np.shape(matrix))
    Px_Cnt_list= []
    for coord in coords:
        row, col = coord
        pxC = Px_Contour(row,col,fam=index_fam[row][col])
        if pxC in Px_Cnt_list:
            for pxCi in Px_Cnt_list:
                if pxC==pxCi:
                    pxC.fam = pxCi.fam
                    break
        else:
            Px_Cnt_list.append(pxC)
        nwpix = pxC.get_bordering_pixels()
        for nwP in nwpix:
            if _acceptable(np.shape(matrix),nwP.coord()):
                if nwP.coord() not in np.array(coords).tolist():
                    nwP.type="bound"
                if nwP in Px_Cnt_list and nwP.type=="mask":
                    for pxCi in Px_Cnt_list:
                        if nwP==pxCi:
                            fam_to_change=pxCi.fam
                            pxCi.fam = nwP.fam
                            break
                    for pxCi in Px_Cnt_list:
                        if pxCi.fam==fam_to_change:
                            pxCi.fam = nwP.fam
                else:
                    Px_Cnt_list.append(nwP)
    return Px_Cnt_list

def extract_family(px_list):
    return list(set([px.fam for px in px_list]))

def interpolate_2D(data,pix,order=0):
    # data:2D matrix
    # pix: coordinate of the pixel to interpolate
    # order: order of interpolation. O=mean between neighboor pixels
    #contours_pix, contours_families_px,contours_families_cont = get_neighboring_coords(data,pix)
    Pix_Cont_list = get_neighboring_coords(data,pix)
    fml_list      =  extract_family(Pix_Cont_list)
    for fam in fml_list:
        vlpx_fm_i  = []
        mskpx_fm_i =[]
        for px in Pix_Cont_list:
            if px.fam==fam:
                row,col = px.coord()
                if px.type=="bound":
                    vlpx_fm_i.append(data[row][col])
                else:
                    mskpx_fm_i.append(px.coord())
        if order==0:
            interp_val = np.nanmean(vlpx_fm_i)
        else:
            #pragma: no cover
            raise
        for i,j in mskpx_fm_i:
            data[i][j] = interp_val
    return data 


def get_numPix(setting,twodim=False):
    setting    = get_setting_module(setting).setting()
    image_path = setting.data_path+setting.image_name
    image      = load_fits(image_path)
    numPix     = image.shape
    if twodim:
        return numPix[0],numPix[1]
    if numPix[0] == numPix[1]:
        return numPix[0]
    else: #pragma: no cover
        raise RuntimeError("Image must be squared, numPix should be indentical in both xy not",numPix)

def get_pixscale(setting_or_tm,only_one=True):
    # check
    # https://danmoser.github.io/notes/gai_fits-imgs.html
    if type(setting_or_tm)!=type(np.array([])) or type(setting_or_tm)!=list:
        transform_pix2angle = get_transf_matrix(setting_or_tm,False)
    else:
        transform_pix2angle = setting_or_tm
    pix_scale1 = np.sqrt(transform_pix2angle[0][0]**2 + transform_pix2angle[1][0]**2)
    pix_scale2 = np.sqrt(transform_pix2angle[0][1]**2 + transform_pix2angle[1][1]**2)
    if only_one:
        if abs(pix_scale1-pix_scale2)>min([pix_scale1,pix_scale2])*.1/100.: # pragma no cover
            # if we are considering only one and the difference is higher then 0.1 % of the smaller pix_scale
            # we raise and error
            raise RuntimeError("Pixel scale significantly different for the x and y axis")
        return pix_scale1
    else: 
        return pix_scale1, pix_scale2
        
def get_rotangle(setting_or_tm,in_deg=True):
    # check
    # https://danmoser.github.io/notes/gai_fits-imgs.html
    # get rotation angle of FoR with respect to pixel coordinates 
    # hence to rotate it correctly we need to rotate the angle by the inverse of this
    if type(setting_or_tm)!=type(np.array([])) and type(setting_or_tm)!=list:
        transform_pix2angle = get_transf_matrix(setting_or_tm,False)
    else:
        transform_pix2angle = setting_or_tm
    rot_angle = np.arctan2(transform_pix2angle[1][0],transform_pix2angle[0][0])
    if not in_deg:
        warnings.warn("Return rotation angle in rad") 
        return rot_angle
    else:
        warnings.warn("Return rotation angle in degrees") 
        return rot_angle*180/np.pi

def extract_phi_ll(model_name,setting,min_ra=0.):
    data_path  = setting.data_path
    trasnforM  = setting.transform_pix2angle
    
    model_path = data_path+"/"+model_name
    param_val  = load_fits(model_path,HDU=-1)

    #param_names = list(param_val.dtype.names)
    pa   = param_val["pa"][1:]
    ra14 = param_val["a14"][1:]
    
    # ignore PA close to the center:
    pa_cut = pa[np.where(np.array(ra14)>min_ra**(1./4))]
    PA     = np.mean(   ) # this should be in the ref. Frame of the image
    # have to be converted first in WST FoR    
    
    rotang    = get_rotangle(trasnforM)
    #phi_lnstr = rotang - PA # see notes 8th June 22
    phi_lnstr = (-PA + rotang) - 180
    return phi_lnstr

def extract_q_ll(model_name,setting,min_ra=0.):
    data_path  = setting.data_path
    trasnforM  = setting.transform_pix2angle
    
    model_path = data_path+"/"+model_name
    param_val  = load_fits(model_path,HDU=-1)

    #param_names = list(param_val.dtype.names)
    eps  = param_val["eps"][1:]
    ra14 = param_val["a14"][1:]
    
    # ignore EPS close to the center:
    eps_cut = eps[np.where(np.array(ra14)>min_ra**(1./4))]
    EPS     = np.mean(eps_cut)
    return EPS


def fits_with_copied_hdr(data,fits_parent_path,data_object="",data_history="",fits_res_namepath=None,overwrite=True,verbose=True):
    
    with fits.open(fits_parent_path,ignore_missing_end=True) as target:
        scihdr  = target[0].header
    
    hdu = fits.PrimaryHDU(data=data,header=scihdr)
    
    if data_object!="":
        hdu.header["OBJECT"]=str(data_object)
    
    if data_history!="":
        hdu.header["HISTORY"]=str(data_history)
    
    if fits_res_namepath is None:
        return hdu
    else:
        if verbose:
            print("saving file "+fits_res_namepath)
        hdu.writeto(fits_res_namepath, overwrite=overwrite)
        return 0
        
def multifits_with_copied_hdr(data_list,fits_parent_path,data_object=[],fits_res_namepath=None,overwrite=True,verbose=True):
    hdu_list = fits.HDUList()
    for i,data in enumerate(data_list):
        if i ==0:
            hdu_list.append(fits_with_copied_hdr(data,fits_parent_path))
            hdr = hdu_list[0].header
        else:
            hdu_list.append(fits.ImageHDU(data=data,header=hdr))
    if data_object != []:
        if len(data_object)!=len(hdu_list):
            raise RuntimeError("Give object name for each data")
        for hdu,dobj in zip(hdu_list,data_object):
            hdu.header["OBJECT"]=str(dobj)
    if fits_res_namepath is None:
        return hdu_list
    else:
        if verbose:
            print("Saving file "+fits_res_namepath)
        hdu_list.writeto(fits_res_namepath, overwrite=overwrite)
        return 0


def _shift_astrometry(data_fits,coord,xy):
    data_fits.header["CRVAL1"]  = coord[0]
    data_fits.header["CRVAL2"]  = coord[1]
    data_fits.header["CRPIX1"]  = xy[0]
    data_fits.header["CRPIX2"]  = xy[1]
    return data_fits

def shift_astrometry(data_fits,coord,xy):
    if type(data_fits) is type(fits.HDUList()):
        shifted_f = fits.HDUList()
        for f in data_fits:
            shifted_f.append(_shift_astrometry(f,coord,xy))
    else:
        shifted_f = _shift_astrometry(data_fits,coord,xy)
    return shifted_f


def get_header(path,hdr_name):
    with fits.open(path,ignore_missing_end=True) as target:
        hdr  = target[0].header
    return float(hdr[hdr_name])

def get_EE(setting,r):
    setting   = get_setting_module(setting).setting()
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
    aper   = [float(aper.split("#")[1]) for aper in titles[2:]]
    diff_r = [aper_i - r for aper_i in aper]
    index_EE = diff_r.index(min(diff_r))
    EE_r = EErow_i[index_EE]
    return float(EE_r)

def get_PHOTFLAM(setting):
    setting   = get_setting_module(setting).setting()
    data_path = setting.data_path
    filt = setting.filter_name.upper()
    if "F1"==filt[:2]:
        # IR
        data_path = "./data/"#"/".join(data_path.split("/")[:-1]) # main dir of all filters
        PF_file = data_path+"/photoflam_IR.csv"
        with open(PF_file,newline="") as f:
                reader = csv.reader(f)
                for i,PFrow in enumerate(reader):
                    if i==0:
                        titles = PFrow
                    if i==1:
                        units  = PFrow
                    if filt in PFrow[0]:
                        PFrow_i = PFrow[1:]
                        break
        index_PF = [titles.index("PHOTFLAM")]
        PF = PFrow_i[index_PF]
        return float(PF)
    else:
        data_path = data_path+"/"+setting.image_name
        return get_header(data_path,"PHOTFLAM")

