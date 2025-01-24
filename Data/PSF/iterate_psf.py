from Utils.tools import *
from Data.image_manipulation import *
from Data.input_data import init_kwrg_data,init_kwrg_numerics,init_kwrg_psf,get_kwargs_model,get_kwargs_constraints
#MOD_PSFITER
from Custom_Model.custom_logL import init_kwrg_custom_likelihood
from Custom_Model.my_lenstronomy.my_fitting_sequence import MyFittingSequence

def _plot_psf_iteration(kw_psf,name,num_iter=None):
    shr     = .4
    psf_in  = kw_psf["kernel_point_source_init"]
    psf_fin = kw_psf["kernel_point_source"]
    #psf_err = kw_psf["psf_error_map"]
    fig,ax  = plt.subplots(1,3)
    ai = ax[0].imshow(np.log10(psf_in.T), cmap = plt.cm.gist_heat)
    fig.colorbar(ai,ax=ax[0],shrink=shr)
    ax[0].set_title("PSF IN")
    ai = ax[1].imshow(np.log10(psf_fin.T), cmap = plt.cm.gist_heat)
    fig.colorbar(ai,ax=ax[1],shrink=shr)
    ax[1].set_title("PSF OUT")
    ai = ax[2].imshow(np.log10(np.abs(psf_fin.T-psf_in.T)), cmap = plt.cm.gist_heat)
    fig.colorbar(ai,ax=ax[2],shrink=shr)
    ax[2].set_title("RESID")
    for a in ax:
        a.axis("off")
    if num_iter:
        fig.suptitle("N* iter.: "+str(num_iter))
    plt.tight_layout()
    plt.savefig(name)

@check_setting
def iterate_psf(setting,kwargs_psf=None,saveplots=False,**kwargs_init_psf):
    if kwargs_psf is None:
        init_kwrg_psf(setting,saveplots=saveplots,**kwargs_init_psf)
    kwargs_data,mask   = init_kwrg_data(setting,saveplots=False,return_mask=True)
    kwargs_numerics    = init_kwrg_numerics(setting)
    kwargs_model       = get_kwargs_model(setting)
    multi_band_list    = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    kwargs_data_joint  = {'multi_band_list': multi_band_list, 
                        'multi_band_type': 'multi-linear'}
    kwargs_constraints = get_kwargs_constraints(setting)
    kwargs_likelihood  = init_kwrg_custom_likelihood(setting,mask,custom=setting.lens_prior.custom_type)
    kwargs_params      = {'lens_model'        : setting.lens_params,
                          'point_source_model': setting.ps_params,
                          'lens_light_model'  : setting.lens_light_params}
    if not setting.WS:
        kwargs_params['source_model'] =   setting.source_params
    psf_iterfitt = MyFittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    psf_iterfitt.psf_iteration(num_iter=getattr(setting,"N_iter_psf",50))
    kwargs_psf   = psf_iterfitt.multi_band_list[0][-2]
    if saveplots:
        savefig_path = get_savefigpath(setting=setting)
        _plot_psf_iteration(kw_psf=kwargs_psf,name=savefig_path+"/psf_iteration.pdf")
    return kwargs_psf