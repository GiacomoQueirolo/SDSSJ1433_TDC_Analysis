#'index_lens_model_list', 'index_source_light_model_list', 'index_lens_light_model_list', 'index_point_source_model_list'

#These arguments should be lists of length the number of imaging bands available and
# each entry in the list is a list of integers specifying the model components being evaluated for the specific band.
# E.g. there are two bands and you want to different light profiles being modeled. 
# - you define two different light profiles 
# lens_light_model_list = ['SERSIC', 'SERSIC'] 
# - set index_lens_light_model_list = [[0], [1]] 
# - (optional) for now all the parameters between the two light profiles are independent in the model. 
# You have the possibility to join a subset of model parameters (e.g. joint centroid). See the Param() class for documentation.

# my note: this means that the different light profile are defined completely independently 
# (eg: main lens light as [SERSIC_ELLIPSE,SERSIC_ELLIPSE] for 2 infrared filters)
# then we defined in each index_..._model_list which band is modelled with what -> f105w ; [0], f140w ;  [1]
# later we then join the required parameters, say the center and ellipticity, with joint_...

from Utils.tools import check_setting,check_if_SUB,check_if_WS

#@check_setting (not needed for this one)
"""
# They should be shared by default
def get_index_lens_model_list(settings):
    # easy, for this is the same for all of them and it's the 3 lens profiles
    index_lens_model =[]
    for _ in settings:
        # Main lens, perturber, shear
        index_lens_model.append([0,1,2])
    return index_lens_model
"""
@check_setting
def get_index_source_light_model_list(settings):
    # this depends on whether there is or there is not a source
    index_source_light_model = []
    ind = 0
    for sett in settings:
        if check_if_WS(sett):
            index_source_light_model.append([])
            #index_source_light_model.append([None]) -> wrong
        else:
            index_source_light_model.append([ind])
            ind+=1
    return index_source_light_model

@check_setting
def get_index_lens_light_model_list(settings):
    # this depends on whether there main lens is subtracted or not
    index_lens_light_model = []
    for sett in settings:
        if check_if_SUB(sett):
            # The perturber and the background are still modelled
            index_lens_light_model.append([1,2])
        else:
            index_lens_light_model.append([0,1,2])
    return index_lens_light_model
    
    
#@check_setting (not needed for this one)
def get_index_point_source_model_list(settings):
    # easy, for this is the same for all of them and it's 1 PS profile
    index_ps_model =[]
    for _ in settings:
        # one PS profile 
        index_ps_model.append([0])
    return index_ps_model

#@check_setting (not needed for this one)
def get_index_point_source_list(settings):
    #point_source_frame_list: list of lists mirroring the structure of the image positions.
    # Integers correspond to the i'th list entry of index_lens_model_list indicating in which frame/band the image is
    # appearing
    # in this case, easy: all band have images
    index_ps_frame_list = []
    for _ in settings:
        index_ps_frame_list.append([0])
    return index_ps_frame_list
