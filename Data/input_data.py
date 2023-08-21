#!/usr/bin/env python
# coding: utf-8

# Define functions to get kwargs_data, kwargs_numerics and kwargs_psf in order to obtain multi_band_list prior/indep.to the model itself

from Utils.tools import *
from Data.image_manipulation import *
from Data.interpolation import interpolate_2D
from lenstronomy.LensModel.lens_model import LensModel
#from Custom_Model.custom_logL import logL_ellipticity_qphi as  logL_ellipticity

@check_setting
def init_kwrg_data(setting,saveplots=False,backup_path="backup_results",return_mask=False):
    if saveplots:
        savefig_path  = get_savefigpath(setting,backup_path)
    image_file   = setting.data_path+setting.image_name
    err_file     = setting.data_path+setting.err_name

    image = load_fits(image_file)
    if saveplots:
        plot_image(image,setting,savefig_path+"/original_image.png")
    err_image = load_fits(err_file)
    # the masking procedure define the mask, subtract the eventual lens light profile,
    # plot the error and image masked and return the subtracted (unmasked) image and the mask
    mask = create_mask(image,setting)
    # Subtract lens light (if necessary)
    image_sub = subtract_light(image,setting)

    if saveplots:
        # Plot the subtracted and masked data image and error
        if setting.sub:
            plot_image(image_sub*mask,setting,savefig_path+"/image_sub_mask.png")
        else:
            plot_image(image_sub*mask,setting,savefig_path+"/image_mask.png")
        plot_image(err_image*mask,setting,savefig_path+"/err_mask.png",err_image=True)

    kwargs_data = { 'noise_map':   err_image,  # noise map
                'ra_at_xy_0':  setting.ra_at_xy_0,  # RA  at (0,0) pixel
                'dec_at_xy_0': setting.dec_at_xy_0, # DEC at (0,0) pixel 
                'ra_shift' :   0., #shifts the coordinate system with respect to 'ra_at_xy_0'    
                'dec_shift':   0., #shifts the coordinate system with respect to 'dec_at_xy_0'    
                'transform_pix2angle': setting.transform_pix2angle,#Mpix2coord,  
                # matrix to translate shift in pixel in shift in relative RA/DEC (2x2 matrix). 
                # Make sure it's units are arcseconds or the angular units you want to model.
                'image_data': image_sub  # 2d data vector
              }
    if not return_mask:
        return kwargs_data
    else:
        return kwargs_data,mask
    
@check_setting
def init_kwrg_psf(setting,saveplots=False,saveinterp=True,backup_path="backup_results"):
    psf_file     = setting.data_path+setting.psf_name 
    err_psf_file = setting.data_path+setting.err_psf
    psf_image = np.array(load_fits(psf_file))
    mask_zero = np.vstack(list(np.where(psf_image<0))).T
    if mask_zero.shape[0]!=0:
        psf_file_interp = psf_file.replace(".fits","_interp.fits") 
        try:
            psf_image = np.array(load_fits(psf_file_interp))
        except FileNotFoundError:
            print("No interpolated PSF found, doing it now")
            psf_image = interpolate_2D(psf_image,mask_zero)
            if saveinterp:
                fits_with_copied_hdr(psf_image,fits_parent_path=psf_file,data_object="PSF image (neg. values interpolated)",data_history="Used image_manipulation.interpolate_2D to interpolate negative values",
                fits_res_namepath=psf_file_interp,overwrite=True,verbose=True)
    
    if saveplots:
        savefig_path  = get_savefigpath(setting,backup_path) 
        if setting.pssf>1:
            plot_image(psf_image,setting,savefig_path+"/psf_supersampled.png")
        else:
            plot_image(psf_image,setting,savefig_path+"/psf.png")
        plot_projection(psf_image,savefig_path)
    #We import the psf error image 
    err_psf_image = load_fits(err_psf_file)
    err_psf_image = psf_correction(err_psf_image,setting)
    if saveplots:
        plot_image(err_psf_image,setting,savefig_path+"/err_psf.png",err_image=True)

    kwargs_psf = {'psf_type': "PIXEL", 
                  'kernel_point_source':psf_image,
                  'point_source_supersampling_factor':setting.pssf,
                  'psf_error_map':err_psf_image}
    return kwargs_psf

@check_setting
def init_kwrg_numerics(setting):
    if setting.pssf==1:
        kwargs_numerics = {'supersampling_factor': 1, 'supersampling_convolution': False}
    else:
        kwargs_numerics = {'supersampling_factor': setting.pssf, 'supersampling_convolution': True,
                            "supersampling_kernel_size":setting.ssks,# int (odd number), size (in regular pixel units) of the super-sampled convolution
                            'point_source_supersampling_factor':setting.pssf} #not sure this is needed -> it is
    return kwargs_numerics

@check_setting
def init_kwrg_likelihood(setting,mask=None):
    kwargs_likelihood = {}
    kwargs_likelihood["check_matched_source_position"]=True
    kwargs_likelihood["source_position_tolerance"]=0.01
    kwargs_likelihood["force_no_add_image"] = True 
    kwargs_likelihood["check_positive_flux"] = True  
    """
    # MOD_CUSTOM_LIKE
    phi_ll = setting.phi_ll if setting.sub else None
    q_ll   = setting.q_ll   if setting.sub else None
    kwargs_likelihood["custom_logL_addition"] = logL_ellipticity(SUB=setting.sub,phi_ll=phi_ll,q_ll=q_ll)
    """
    if mask is None:
        _,mask=init_kwrg_data(setting,saveplots=False,return_mask=True)
    if str(type(mask))=="<class 'numpy.ndarray'>":
        mask=mask.tolist()
    kwargs_likelihood["image_likelihood_mask_list"] =  [mask]
    return kwargs_likelihood

@check_setting
def init_multi_band_list(setting,saveplots=False,backup_path="backup_results",return_mask=False):
    # setting can be the string with the name of the setting or the instance of the relative setting function 
    kwargs_psf      = init_kwrg_psf(setting,saveplots,backup_path)
    kwargs_numerics = init_kwrg_numerics(setting)
    if return_mask:
        kwargs_data,mask = init_kwrg_data(setting,saveplots,backup_path,return_mask=return_mask)
        multi_band_list  = [[kwargs_data, kwargs_psf, kwargs_numerics]]
        return multi_band_list,mask
    else:
        kwargs_data     = init_kwrg_data(setting,saveplots,backup_path,return_mask=False)
        multi_band_list = [[kwargs_data, kwargs_psf, kwargs_numerics]]
        return multi_band_list

@check_setting
def init_lens_model_list(setting):
    lens_model_list = ['SIE']
    if setting.CP:
        #lens_model_list = ['PEMD']
        lens_model_list = ['EPL_NUMBA'] #test
    lens_model_list = [*lens_model_list,'SIS','SHEAR_GAMMA_PSI']
    return lens_model_list

@check_setting
def init_lens_light_model_list(setting):
    if setting.sub==False:
        light_model_list = ["SERSIC_ELLIPSE", "SERSIC","UNIFORM"]
    else:
        light_model_list = ["SERSIC","UNIFORM"]
    if hasattr(setting,"no_pert"):
        if setting.no_pert==True:
            light_model_list=["UNIFORM"]
    if getattr(setting,"no_bckg",False):
        light_model_list.remove("UNIFORM")
    return light_model_list

@check_setting
def init_source_model_list(setting):
    if not setting.WS:
        if getattr(setting,"ellipt",False):
            source_model_list = ['SERSIC_ELLIPSE']
        elif getattr(setting,"shapelets",False):
            source_model_list = ['SHAPELETS']
        # MOD_2SRC: add a second source Sersic
        elif getattr(setting,"second_source",False):
            source_model_list = ['SERSIC','SERSIC']
        else:
            source_model_list = ['SERSIC']
    else:
        source_model_list = None
    return source_model_list


@check_setting
def get_kwargs_model(setting):
    lens_model_list   = init_lens_model_list(setting)
    # Lens light profile with perturber
    light_model_list  =  init_lens_light_model_list(setting)
    # Source host gal
    source_model_list = init_source_model_list(setting)
    
    # Source light model: point 
    point_source_list = ['LENSED_POSITION']

    kwargs_model = {'z_lens':setting.z_lens,
                'z_source':setting.z_source,
                'lens_model_list': lens_model_list,
                'lens_light_model_list': light_model_list,
                'point_source_model_list': point_source_list,
                'additional_images_list': [False], #list of bools (same length as point_source_type_list).
                # If True, search for additional images of the same source is conducted.
                 'fixed_magnification_list': setting.fixed_mag,  # list of bools (same length as point_source_type_list).
                #If True, magnification ratio of point sources is fixed to the one given by the lens model 
                }
    if not setting.WS:
        kwargs_model['source_light_model_list'] = source_model_list
    return kwargs_model

@check_setting
def get_kwargs_constraints(setting):
    if setting.sub ==False:
        # MOD_XYLL and MOD_PLL
        if getattr(setting,"xyll",False) or getattr(setting,"pll",False):
            joint_lens_with_light=[[1,1,["center_x","center_y"]]]
            #  'joint_lens_with_light': list [[i_light, k_lens, ['param_name1', 'param_name2', ...]], [...], ...],
            #   joint parameter between lens model and lens light model
        else:
            joint_lens_with_light=[[0,0,["center_x","center_y"]],[1,1,["center_x","center_y"]]]
    else:
            joint_lens_with_light=[[0,1,["center_x","center_y"]]]
            
    kwargs_constraints = {'num_point_source_list': [4], 
                          'solver_type': 'NONE',
                          'joint_lens_with_light':joint_lens_with_light}
    # mod free source
    if not check_if_WS(setting) and not setting.FS:
        kwargs_constraints['joint_source_with_point_source'] = [[0, 0]]
    return kwargs_constraints

@check_setting
def init_lens_model(setting,lens_model_list=None):
    if lens_model_list is None:
        lens_model_list = init_lens_model_list(setting)
    lensModel = LensModel(lens_model_list=lens_model_list, z_lens=setting.z_lens, z_source=setting.z_source)
    return lensModel