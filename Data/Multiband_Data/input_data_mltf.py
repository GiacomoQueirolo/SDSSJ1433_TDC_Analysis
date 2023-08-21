#!/usr/bin/env python
# coding: utf-8

# Define functions to get kwargs_data, kwargs_numerics and kwargs_psf in order to obtain multi_band_list prior/indep.to the model itself

from Utils.tools import *
from Utils.Multiband_Utils.index_models import *
from Utils.Multiband_Utils.tools_multifilter import *
from Utils.Multiband_Utils.tools_multifilter import _append
from Utils.Multiband_Utils.Setting_multifilter import *

from Data.image_manipulation import *
from Data.input_data import init_kwrg_data,init_kwrg_psf,init_kwrg_numerics,init_lens_model_list

from lenstronomy.LensModel.lens_model import LensModel

@check_mltf_setting
def init_lens_model_list_mltf(multifilter_sett):
    # these stay the same no matter how many filters we consider
    # check index_models.py initial comment 
    # Test that the single setting files are compatible
    lens_model_list = init_lens_model_list(multifilter_sett.settings[0])
    for sett_i in multifilter_sett.settings[1:]:
        if lens_model_list!= init_lens_model_list(sett_i):
            raise RuntimeError(f"{get_setting_name(sett_i)} has different lens model list than the others")
    """lens_model_list = ['SIE']
    if multifilter_sett.CP:
        #lens_model_list = ['PEMD']
        lens_model_list = ['EPL_NUMBA'] #test
    lens_model_list = [*lens_model_list,'SIS','SHEAR_GAMMA_PSI']"""
    return lens_model_list

@check_mltf_setting
def init_lens_light_model_list_mltf(multifilter_sett):
    # for each filter we need to define a light model
    light_model_list = []
    for sett in multifilter_sett.settings:
        if check_if_SUB(sett):
            _append(light_model_list,["SERSIC","UNIFORM"])
        else:
            _append(light_model_list,["SERSIC_ELLIPSE", "SERSIC","UNIFORM"])
    return light_model_list



@check_mltf_setting
def init_source_model_list_mltf(multifilter_sett):
    source_model_list = []
    for sett in multifilter_sett.settings:
        if not check_if_WS(sett):
            if hasattr(sett,"ellipt"):
                if sett.ellipt:
                    source_model_list.append('SERSIC_ELLIPSE')
                else:
                    source_model_list.append('SERSIC')
            #test w. shapelts
            elif hasattr(sett,"shapelets"):
                if sett.shapelets:
                    source_model_list.append('SHAPELETS')
                else:
                    source_model_list.append('SERSIC')
            else:
                source_model_list.append('SERSIC_ELLIPSE')
    return source_model_list


@check_mltf_setting
def get_kwargs_model_mltf(multifilter_sett):
    lens_model_list   = init_lens_model_list_mltf(multifilter_sett)
    # Lens light profile with perturber
    light_model_list  =  init_lens_light_model_list_mltf(multifilter_sett)
    # Source host gal
    source_model_list = init_source_model_list_mltf(multifilter_sett)
    
    # Source light model: point 
    point_source_list = ['LENSED_POSITION']

    kwargs_model = {'z_lens':multifilter_sett.z_lens,
                    'z_source':multifilter_sett.z_source,
                    'lens_model_list': lens_model_list,
                    'lens_light_model_list': light_model_list,
                    'point_source_model_list': point_source_list,
                    'additional_images_list': [False], #list of bools (same length as point_source_list).
                    # If True, search for additional images of the same source is conducted.
                    'fixed_magnification_list': multifilter_sett.fixed_mag,  # list of bools (same length as point_source_list).
                    #If True, magnification ratio of point sources is fixed to the one given by the lens model 
                    }
    if not multifilter_sett.allWS:
        kwargs_model['source_light_model_list'] = source_model_list
    
    #kwargs_model["index_lens_model_list"]         = get_index_lens_model_list(multifilter_sett.settings)
    kwargs_model["index_lens_light_model_list"]   = get_index_lens_light_model_list(multifilter_sett.settings)
    kwargs_model["index_source_light_model_list"] = get_index_source_light_model_list(multifilter_sett.settings)
    kwargs_model["index_point_source_model_list"] = get_index_point_source_model_list(multifilter_sett.settings)
    #point_source_frame_list: list of lists mirroring the structure of the image positions.
    # Integers correspond to the i'th list entry of index_lens_model_list indicating in which frame/band the image is
    # appearing
    kwargs_model['point_source_frame_list']       = get_index_point_source_frame_list(multifilter_sett.settings)
    return kwargs_model

@check_mltf_setting
def init_kwrg_likelihood_mltf(multifilter_sett,masks=None):
    kwargs_likelihood = {}
    kwargs_likelihood["check_matched_source_position"]=True
    kwargs_likelihood["source_position_tolerance"]=0.01
    kwargs_likelihood["force_no_add_image"] = True 
    kwargs_likelihood["check_positive_flux"] = True  
    kwargs_likelihood["image_likelihood_mask_list"] = []
    for ist,sett in enumerate(multifilter_sett.settings):
        if masks is None:
            _,mask=init_kwrg_data(sett,saveplots=False,return_mask=True)
        else:
            mask = masks[ist]
        if str(type(mask))=="<class 'numpy.ndarray'>":
            mask=mask.tolist()
        kwargs_likelihood["image_likelihood_mask_list"].append(mask)
    return kwargs_likelihood

@check_mltf_setting
def get_kwargs_constraints_mltf(multifilter_sett,kwargs_model):
    kwargs_constraints = {'solver_type': 'NONE',
                          'num_point_source_list':[],
                          'joint_lens_with_light':[],
                          }
    if len(multifilter_sett.settings)!=1:
        kwargs_constraints['joint_lens_with_lens']=[]
    if not multifilter_sett.allWS:
        kwargs_constraints['joint_source_with_point_source']=[]
    #source_list        = kwargs_model["source_light_model_list"]
    lens_light_list     = kwargs_model["lens_light_model_list"]
    indexes_pert_light  = np.where(np.array(lens_light_list)=="SERSIC")[0]
    indexes_mainL_light = np.where(np.array(lens_light_list)=="SERSIC_ELLIPSE")[0]

    mll_ind,pl_ind = 0, 0 
    ind_sources    = 0
    for setting in multifilter_sett.settings:
        kwargs_constraints['num_point_source_list'].append(4)
        if not check_if_SUB(setting):
            # MOD_XYLL and MOD_PLL
            if getattr(setting,"xyll",False) or getattr(setting,"pll",False):
                joint_lens_with_light=[indexes_mainL_light[mll_ind],1,["center_x","center_y"]]
                mll_ind +=1
                joint_lens_with_light=[indexes_pert_light[pl_ind],1,["center_x","center_y"]]
                pl_ind+=1
                #  'joint_lens_with_light': list [[i_light, k_lens, ['param_name1', 'param_name2', ...]], [...], ...],
                #   joint parameter between lens model and lens light model
            else:
                raise RuntimeWarning("Pragma no cover: not yet implemented for not pll settings")
                #joint_lens_with_light=[[0,0,["center_x","center_y"]],[1,1,["center_x","center_y"]]]
        else:
            joint_lens_with_light=[indexes_pert_light[pl_ind],1,["center_x","center_y"]]
            pl_ind+=1

        if not check_if_WS(setting) and not setting.FS:
            #joint_source_with_point_source : list [[i_point_source, k_source], [...],
            kwargs_constraints['joint_source_with_point_source'].append([0, ind_sources])
            ind_sources+=1

        kwargs_constraints['joint_lens_with_light'].append(joint_lens_with_light)
    # now we should define the constraints BETWEEN the filters:
    # point sources are modelled with the same profile, so no need to constrain it differently
    # > is it a problem due the frame of reference?
    # sources -> fixed to the PS, other parameter should be independent
    # lenses > all same profiles
    # lensed light -> pos of perturber is constrained from lens mass -> done
    # uniform bck is independent
    # main lens light: -> center and ellipticities might be fixed?
    # ellipticities: no, they might change over the colors
    # center: yes-> [[i_ml_light,j_ml_light,["center_x","center_y"]]]       
    
    if not multifilter_sett.allSUB:
        kwargs_constraints['joint_lens_light_with_lens_light']=[]
        # we want to link each of them 
        # all to the first model
        first_ll_index =  lens_light_list.index("SERSIC_ELLIPSE")
        for i,llp in enumerate(lens_light_list):
            if i!=first_ll_index and llp=="SERSIC_ELLIPSE":
                kwargs_constraints['joint_lens_light_with_lens_light'].append([first_ll_index,i,["center_x","center_y"]]) 
                    
    return kwargs_constraints

@check_mltf_setting
def get_fixed_sources_list_mltf(multifilter_sett):
    # we sett all params of the sources to be fixed
    # for the first iteration of the model
    fixed_sources_list = []
    nsources = 0
    for sett in multifilter_sett.settings:
        if not check_if_WS(sett):
            # all sources should be Sersic
            # the center cannot be fixed bc that is connnected to the images
            fixed_sources_list.append([nsources,['R_sersic', 'n_sersic']])
    return fixed_sources_list

@check_mltf_setting
def init_multi_band_list_mltf(multifilter_sett,saveplots=False, return_masks=False):
    # setting can be the string with the name of the setting or the instance of the relative setting function 
    multi_band_list = []
    masks = []
    for sett in multifilter_sett.settings:
        kwargs_psf      = init_kwrg_psf(sett,saveplots,multifilter_sett.backup_path)
        kwargs_numerics = init_kwrg_numerics(sett)
        if return_masks:
            kwargs_data,mask = init_kwrg_data(sett,saveplots,multifilter_sett.backup_path,return_mask=return_masks)
            multi_band_list.append([kwargs_data, kwargs_psf, kwargs_numerics])
            masks.append(mask)
        else:
            kwargs_data     = init_kwrg_data(sett,saveplots,multifilter_sett.backup_path,return_mask=False)
            multi_band_list.append([kwargs_data, kwargs_psf, kwargs_numerics])
    if return_masks:
        return multi_band_list,masks
    return multi_band_list



@check_mltf_setting
def init_lens_model(multifilter_sett,lens_model_list=None):
    if lens_model_list is None:
        lens_model_list = init_lens_model_list_mltf(multifilter_sett)
    lensModel = LensModel(lens_model_list=lens_model_list, z_lens=multifilter_sett.z_lens, z_source=multifilter_sett.z_source)
    return lensModel
