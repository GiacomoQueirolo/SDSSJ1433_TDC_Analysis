# Reworked copy of Rewamp_MtL_post; now considering colour, both filters and is independent from previous versions

import sys
import argparse
import numpy as np
from corner import quantile 
from acstools import acszpt 
from astropy import units as u
import matplotlib.pyplot as plt
from astropy import constants as const
from photutils.aperture import EllipticalAperture


from Utils.tools import *
from Utils.get_res import *
from Utils.conversion import *
from MassToLight.grid_class import dist
from Data.input_data import init_lens_model
from Utils.statistical_tools import simp,easy_res_str
from Posterior_analysis.tools_Post import default_cosmo
from Data.image_manipulation import load_fits,get_rotangle

#####################################
# Geom. and cosmo. functions
####################################

def get_annulus(r,eps,from_0=True):
    # r is considered in arcsec or px, not **1/4
    # copied from sersic_def.py of Matthias Kluge
    area_ann = np.pi*(r[1:]**2*(1-eps[1:]) - r[:-1]**2*(1-eps[:-1]) )
    if from_0:
        area_ann = [np.pi*(r[0]**2)*(1-eps[0]), *area_ann ]
    try:
        unit_area = area_ann[0].unit
        return np.array([an.value for an in area_ann])*unit_area
    except AttributeError:
        return np.array(area_ann)


@check_setting
def get_cosmo_prm(setting,cosmo=default_cosmo,SigCr_arcs2=False):
    # SigCr_arcs2: if True, convert Sigma_Crit in M_sun/arcsec^2
    # instead of M_sun/kpc_proper^2
    z_source  = setting.z_source
    z_lens    = setting.z_lens 
    # Angular diameter distance [...] gives the proper (sometimes called ‘physical’) transverse distance 
    # corresponding to an angle of 1 radian for an object at redshift z.
    cosmo_dd  = cosmo.angular_diameter_distance(z_lens).to("kpc")   #kpc
    cosmo_ds  = cosmo.angular_diameter_distance(z_source).to("kpc") #kpc
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source).to("kpc") #kpc
    
    # Sigma_Crit = D_s*c^2/(4piG * D_d*D_ds)
    Sigma_Crit   = (cosmo_ds*const.c**2)/(4*np.pi*const.G*cosmo_dds*cosmo_dd)
    Sigma_Crit   = Sigma_Crit.to("Msun /(kpc kpc)")    
    cosmo_Ld     = cosmo.luminosity_distance(z_lens) 
    if SigCr_arcs2:
        Sigma_Crit/=(cosmo.arcsec_per_kpc_proper(z_lens)**2)
    return {"cosmo":cosmo,"Sigma_Crit":Sigma_Crit,"cosmo_dd":cosmo_dd,"cosmo_Ld":cosmo_Ld,"arcsec_per_kpc_proper":cosmo.arcsec_per_kpc_proper(z_lens)}

######################################
######################################
    
###################
# Mass functions
###################


def get_kappa(r,pa,kwres,lens_model,transform_pix2angle):#,phi_thresh=5,q_thresh=.1):
    # have to be centered on the lens center
    x_0   = kwres["kwargs_lens"][0]["center_x"]*u.arcsec # those are good for both
    y_0   = kwres["kwargs_lens"][0]["center_y"]*u.arcsec # lens and mass models
    # The light is measured based on the annulus, which is in turn defined on the ellipticity of the isoph-light model
    # we have to measure the mass in the same way -> but indeed this is then multiplied to an elliptical annulus
    # which might differ from the PA and axis ratio of the model 
    # A) assume that k at r,eps -> x,y is almost identical to the mean of k around the annulus
    # B) compute the mean -> don't

    # r is assumed to be in arcsecs
    # pa is assumed to be in deg -> to convert in our reference point
    rotang = get_rotangle(transform_pix2angle)*u.degree
    phi    = (rotang-pa) - 180*u.degree
    if np.mean(phi)<-100*u.degree:
        phi+=180*u.degree
    
    x = x_0 + r*np.cos(phi) 
    y = y_0 + r*np.sin(phi) 

    kappa = lens_model.kappa(x=x.to("arcsec").value,y=y.to("arcsec").value,kwargs=kwres['kwargs_lens'])
    return kappa

    
def _get_cumulative_mass(kappa,r,eps,crit_density,verbose=False):
    #if getattr(crit_density,"value",False):
    #    print("Crit density {:.2e}".format(crit_density))
    #    crit_density = crit_density.to("Msun/(arcsec arcsec)").value
    try:
        crit_density.value
    except AttributeError:
        crit_density = crit_density*u.Msun/(u.arcsec**2) 
    area_annulus = get_annulus(r=r,eps=eps,from_0=True) #area in arcsec**2
    dens   = kappa*crit_density
    mass_i = area_annulus*dens
    cumulative_mass = [mass_i[0]]
    for i in range(1,len(kappa)):
        cumulative_mass.append(mass_i[i]+cumulative_mass[-1])
    if verbose:
        print("average kappa ",np.mean(kappa))
        print("total area ",np.sum(area_annulus))
    unit_cm = cumulative_mass[0].unit
    return np.array([cm.value for cm in cumulative_mass])*unit_cm


def get_cumulative_mass(sett,r,eps,pa,kwres=None,lens_model=None,cosmo=default_cosmo,verbose=False):
    if lens_model is None:
        lens_model = init_lens_model(sett)
    if kwres is None:
        kwres = get_kwres(sett)["kwargs_results"]
    kappa           = get_kappa(r,pa,kwres=kwres,lens_model=lens_model,transform_pix2angle=sett.transform_pix2angle)
    Sigma_Crit      = get_cosmo_prm(sett,cosmo=cosmo,SigCr_arcs2=True)["Sigma_Crit"]
    cumulative_mass = _get_cumulative_mass(kappa=kappa,r=r,eps=eps,crit_density=Sigma_Crit,verbose=verbose)
    return cumulative_mass

def M_E(setting,kw_res=None):
    # output the Einstein Mass (avrg. mass at the Einstein radius)
    if kw_res is None:
        kw_res = getattr(sett,"kwres",get_kwres(setting)["kwargs_results"])
    lens_res = kw_res["kwargs_lens"][0]
    theta_E  = lens_res["theta_E"]*u.arcsec.to("rad")
    kw_cosmo      = get_cosmo_prm(setting)
    Sigma_Crit    = kw_cosmo["Sigma_Crit"]
    
    cosmo_dd      = kw_cosmo["cosmo_dd"]
    Mass_Einstein = np.pi*(theta_E**2)*(cosmo_dd**2)*Sigma_Crit
    return Mass_Einstein


def velocity_disp_theta(theta,z_lens,z_source,cosmo=default_cosmo):
    cosmo_ds  = cosmo.angular_diameter_distance(z_source)
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source)
    theta_rad = np.pi*theta/(180*3600)
    return const.c.to("km/s")*np.sqrt(theta_rad*cosmo_ds/(4*np.pi*cosmo_dds))


######################################
######################################
    
#####################
# Light functions
######################


def get_K_corr(setting):
    filter = get_filter(setting)
    # hardcoded from http://kcor.sai.msu.ru/
    # z_lens = 0.407
    # Johnson B, color B-Ic =1.2
    if "475" in filter:
        return 0.09 
    # Cousins Ic , color V-Ic =0.5 
    # WEIRD: to obtain -0.14, we have to give V-Ic=0.7
    elif "814" in filter:
        return -0.14
    else:
        #pragma no cover
        raise RuntimeError("Not implemented")
    
@check_setting
def get_gal_extinction(setting):
    filter = get_filter(setting)
    # hardcoded from 
    # https://ned.ipac.caltech.edu/extinction_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=14%3A33%3A26.5602&dec=%2B60%3A07%3A42.858
    # RA  = 14:33:26.5602
    # DEC = +60:07:42.858
    if "475" in filter:
        return 0.029
    # Cousins Ic , color V-Ic =0.5
    elif "814" in filter:
        return 0.014
    else:
        #pragma no cover
        raise RuntimeError("Not implemented")
    
def get_ZP(filt):
    try:
        if "X" in filt:
            filt = filt.replace("X","W")
        q = acszpt.Query(date='2017-01-01', detector='WFC', filt=filt)
        filter_zpt = q.fetch()
        ZP=filter_zpt['ABmag'][0].value
    except TypeError:
        print("The automatic query failed, using harcoded ZP")
        if "475" in filt:
            ZP = 26.044 # AB_mag
        elif "814" in filt:
            ZP = 25.937 #AB_mag
        else:
            raise RuntimeError(f"No ZP defined for filter {filt}")
    return ZP

def get_M_sun(filt):
    if "475" in filt:
        M_sun=5.07# AB_mag
    elif "814" in filt:
        M_sun=4.52 #AB_mag 
    else:
        raise RuntimeError(f"No M_sun defined for filter {filt}")
    return M_sun

def get_error_flux(x0,y0,a,eps,pa,data,err_frame):
    flux,err_flux = [],[]
    for i in range(len(x0)):
        aper = EllipticalAperture((x0[i],y0[i]),a[i],a[i]*(1-eps[i]),theta=np.deg2rad(pa[i]))
        fl,erflx= aper.do_photometry(data,err_frame)
        flux.append(fl)
        err_flux.append(erflx)
    return flux,err_flux

def get_error_mag(flux,err_flux,pix_scale):
    errmag =[] # AB mag
    errMag =[] # AB mag
    errLum =[] # solar_lum
    
    for i in range(len(err_flux)):
        flx_pix     = flux[i]/(pix_scale**2)
        err_flx_pix = err_flux[i]/(pix_scale**2)
        errLum.append(err_flx_pix[0])
        errm        = np.std(-2.5*np.log10(np.random.normal(flx_pix,err_flx_pix,10000)))
        #errm = 2.5*err_flx_pix/(flx_pix*np.log(10))
        errmag.append(errm)
        errMag.append(errm)
    errmag = np.array(errmag)*u.mag
    errMag = np.array(errMag)*u.mag
    errLum = np.array(errLum)*u.solLum
    return errLum,errMag,errmag

def get_luminosity(flux,M_sun,distance,ZP,pix_scale,z_lens,K,gal_ext):
    mag =[] # AB mag
    Mag =[] # AB mag
    Lum =[] # solar_lum

    d_mod=5*np.log10(distance.to("pc")/(10*u.pc))
    d_mod=d_mod.value
    cosmic_dimming = 2.5*np.log10((1+z_lens)**4)
    
    for i in range(len(flux)):
        flx_pix = flux[i]/(pix_scale**2)
        try:
            mag.append(-2.5*np.log10(flx_pix.value)+ZP)
        except AttributeError:
            mag.append(-2.5*np.log10(flx_pix)+ZP)
        Mag.append(mag[i]-d_mod-K-cosmic_dimming - gal_ext)#consider galactic extinction , cosmic dimming, Kcorr
        Lum.append(10**(0.4*(M_sun-Mag[i])))
    mag = np.array(mag)*u.mag
    Mag = np.array(Mag)*u.mag
    Lum = np.array(Lum)*u.solLum
    return Lum,Mag,mag

def new_get_luminosity(sett,cosmo=default_cosmo):
    model_name = sett.lens_light_model_name_full
        
    filt = get_filter(sett)
    if "475" in filt:
        # problem: err_file is cut wr to the full image
        # hence we have to correct the center coords
        cut_x,cut_y = (180,336),(160,316)
        i_hdu = 8
    elif "814" in filt:
        model_name = "isopy_220527.fits"
        cut_x,cut_y = (833,989),(844,1000)
        i_hdu = 4
    else:
        raise RuntimeWarning("Filter not recognised:"+filt)
    image_path = sett.data_path+"/"+model_name 
    try:
        sett.rad
    except:
        a14 = load_fits(image_path,HDU=-1).field("a14")[1:]
        sett.rad  = (a14**4)*u.arcsec # this is actually a, ie the major axis
    a           = sett.rad/sett.pix_scale
    data        = load_fits(image_path,HDU=i_hdu)[cut_y[0]:cut_y[1]+1,cut_x[0]:cut_x[1]+1]
    error_frame = load_fits(sett.data_path+sett.err_name)      
    flux , err_flux = get_error_flux(sett.x0,sett.y0,a=a.value,
                eps=sett.eps_light,pa=sett.pa_light.value,data=data,
                err_frame=error_frame)  
    K                    = get_K_corr(setting=sett)
    ZP                   = get_ZP(filt=filt)
    M_sun                = get_M_sun(filt=filt)
    gal_ext              = get_gal_extinction(sett)
    cosmo_ld             = cosmo.luminosity_distance(z=sett.z_lens).to("pc")
    Lum,Mag,mag          = get_luminosity(flux,M_sun,cosmo_ld,ZP,sett.pix_scale.value,sett.z_lens,K,gal_ext)
    errLum,errMag,errmag = get_error_mag(flux,err_flux,sett.pix_scale.value)
    return Lum,Mag,mag,errLum,errMag,errmag


def get_effective_rad(rad,lum_sum):
    # simply the radius at which the luminosity is half the total luminosity
    arr_min      = np.abs(lum_sum-(np.max(lum_sum)/2))
    ind_half_lum = int(np.where(arr_min==np.min(arr_min))[0])
    return rad[ind_half_lum]

def produce_color_plots(sett1,sett2,max_rad=3):
    new_rad  = sett1.rad.value[np.where(sett1.rad.value<max_rad)]
    mag1 = sett1.mag.value[:len(new_rad)]
    mag2 = sett2.mag.value[:len(new_rad)]
    errmag1,errmag2 = sett1.errmag.value[:len(new_rad)],sett2.errmag.value[:len(new_rad)]
    
    for st,mg,errmg in [(sett1,mag1,errmag1),(sett2,mag2,errmag2)]:
        plt.scatter (new_rad[2:],mg[2:])
        plt.errorbar(new_rad[2:],mg[2:],yerr=errmg[2:],fmt="b+")
        plt.savefig(f"{get_savefigpath(st)}/new_mag.pdf")
        plt.close()
    
    col     = mag1 - mag2
    err_col = np.array([np.sqrt((em1**2)+(em2**2)) for em1,em2 in zip(errmag1,errmag2)])
 
    r,d   = [sett1.ps_params[0][0].get(key+"_image") for key in ["ra","dec"]]
    cx,cy = [sett1.lens_prior().pll[k][0] for k in ["x","y"]]
    dists   = [dist([ri,di],[cx,cy]) for ri,di in zip(r,d)]
    im0,im1 = min(dists),max(dists)
    regim_i = np.where((new_rad[2:]>im0) & (new_rad[2:]<im1))
    mean_col_im = np.array(col)[regim_i].mean()
    print("Mean color at images loc ",np.round(mean_col_im,2))
    print("Color at rad ",new_rad[-1],"=",col[-1],"\pm",err_col[-1])
    plt.axvspan(0,0.05,alpha=.5,color="grey",label="PSF Dominated Area")
    plt.errorbar(new_rad[3:],col[3:],yerr=err_col[3:],fmt="b.")
    plt.plot(new_rad[3:],col[3:],".")
    plt.xlabel(r'a ["]')
    plt.ylabel(r"m$_{"+get_filter(sett1).upper()+r"}$-m$_{"+get_filter(sett2).upper()+r"}$ [mag]")
    plt.axvspan(im0,im1,alpha=.5,color="green",label="Area of images location")
 
    
    plt.title("Main Lens Colour Profile")
    plt.legend()
    for sett in setts:
        color_plot_name = f"{get_savefigpath(sett)}/color_{get_filter(sett1).upper()}_{get_filter(sett2).upper()}_interpolated.pdf"
        plt.savefig(color_plot_name)
        print(f"Saved {color_plot_name}") 
    plt.close() 
######################################
######################################
    
#######################
# Resulting  functions
#######################

def get_M_w_R(sett,mcmc=None,param_mcmc=None,cosmo=default_cosmo,ret_w_units=False,cut_mcmc=False):
    m_cumul        = []
    M_E_mcmc       = []
    M_theta_mcmc   = []
    theta_E_mcmc   = []
    sig_theta_mcmc = []

    if mcmc is None:
        mcmc = get_mcmc_smpl(sett)
    if param_mcmc is None:
        param_mcmc = get_mcmc_prm(sett)
    if cut_mcmc:
        print("initial length mcmc ",len(mcmc))
        fact_cut  = 1/2
        fact_smpl = 1/20
        n_cut  = int(len(mcmc)*fact_cut)
        mcmc   = mcmc[n_cut:]
        n_smpl = int(len(mcmc)*fact_smpl)
        mcmc   = np.array(mcmc)[np.random.choice(np.arange(len(mcmc)),n_smpl)]
        print("cut length mcmc     ",len(mcmc))
    
    for i in range(len(mcmc)):
        if i%(int(len(mcmc)/10))==0:
            print(np.round(i*100/len(mcmc),0),"%")
        kwargs_result_i    = sett.produce_kwargs_result(mcmc,param_mcmc,i)
        sigma_theta        = velocity_disp_theta(kwargs_result_i["kwargs_lens"][0]["theta_E"],sett.z_lens,sett.z_source,cosmo)
        M_Ei               = M_E(sett,kw_res=kwargs_result_i)
        q_mass_i,pa_mass_i = qphi_from_e1e2(e1=kwargs_result_i["kwargs_lens"][0]["e1"],e2=kwargs_result_i["kwargs_lens"][0]["e2"])
        pa_mass_i   *= u.degree
        eps_mass_i   = np.ones_like(sett.rad.value)*(1-q_mass_i)
        mass_cumul_i = get_cumulative_mass(sett=sett,r=sett.rad,eps=eps_mass_i,pa=pa_mass_i,kwres=kwargs_result_i,\
                                            lens_model=lens_model,cosmo=cosmo)
        theta_E_i = kwargs_result_i["kwargs_lens"][0]["theta_E"]*u.arcsec
        i_theta_E = int(np.where(np.abs(theta_E_i-sett.rad)==np.min(np.abs(theta_E_i-sett.rad)))[0])
        M_theta_i = mass_cumul_i[i_theta_E]
        if i==0:
            unit_M_theta_mcmc   = M_theta_i.unit
            unit_M_E_mcmc       = M_Ei.unit
            unit_sig_theta_mcmc = sigma_theta.unit
        M_E_mcmc.append(M_Ei.value)
        m_cumul.append(mass_cumul_i.value)
        M_theta_mcmc.append(M_theta_i.value)
        theta_E_mcmc.append(theta_E_i.value)
        sig_theta_mcmc.append(sigma_theta.value)
    if ret_w_units:
        return np.array(sig_theta_mcmc)*unit_sig_theta_mcmc,np.array(M_E_mcmc)*unit_M_E_mcmc,\
            np.array(M_theta_mcmc)*unit_M_theta_mcmc,np.array(m_cumul)*unit_M_theta_mcmc,np.array(theta_E_mcmc)*u.arcsec
    else:
        sig_theta_mcmc,M_E_mcmc,M_theta_mcmc,m_cumul,theta_E_mcmc



    

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(prog=sys.argv[0],description="Mass with posterior computed (note: tailored for J1433 analysis, only viable for optical filters)",formatter_class=CustomFormatter)
    parser.add_argument("--check_prior_res", action="store_true", dest="check_prior_res", default=False,
                    help="Check if prior results are available, if so loads them")
    parser.add_argument("--cut_mcmc", action="store_true", dest="cut_mcmc", default=False,
                    help="Cut mcmc to a shorter version")
    args      = parser.parse_args()
    cut_mcmc  = args.cut_mcmc
    prior_res = args.check_prior_res
    cosmo     = default_cosmo
    sett_f475 = get_setting_module("f475w_CPRP_PLLI_ws",1)
    sett_f814 = get_setting_module("f814w_CPRP_PLLI_ws",1)
    setts     = [sett_f475,sett_f814]


    fnt = 14
    plt.rcParams['xtick.labelsize'] = fnt
    plt.rcParams['ytick.labelsize'] = fnt 
    plt.rcParams['font.size'] = fnt
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('font', size=fnt)          # controls default text sizes
    plt.rc('axes', titlesize=fnt)     # fontsize of the axes title
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('legend', fontsize=fnt)    # legend fontsize


    ######################################
    # Mass and MtL ratio
    lens_model = init_lens_model(setting=sett_f814) #should be the same for both sett

    for sett in setts: 
        svpth = get_savemcmcpath(sett)
        # to better define
        sett.kwres = get_kwres(sett)["kwargs_results"]
        data_path  = sett.data_path
        model_name = sett.lens_light_model_name_full
        if "814" in get_filter(sett):
            model_name = "isopy_220527.fits"
            # problem: err_file is cut wr to the full image
            # hence we have to correct the center coords
            cut_x,cut_y = 833,844 
        else:
            cut_x,cut_y = 180,160 
        model_path = f"{data_path}/{model_name}"

        param_val = load_fits(model_path,HDU=-1)
        a14       = param_val.field("a14")[1:]
        sett.rad  = (a14**4)*u.arcsec # this is actually a, ie the major axis          
        sett.pix_scale = sett.pix_scale*u.arcsec/u.pixel 
        sett.x0   = np.array(param_val.field("x0")[1:])- cut_x 
        sett.y0   = np.array(param_val.field("y0")[1:])- cut_y
        sett.eps_light = param_val.field("eps")[1:]
        sett.pa_light  = param_val.field("pa")[1:]*u.degree
        
        #######################
        sett.Lum,sett.Mag,sett.mag,sett.errLum,sett.errMag,sett.errmag = new_get_luminosity(sett)
        sett.Lum = sett.Lum[:,0]
        ########################
        svnm_sig_theta = save_json_name(setting=sett,path=svpth,filename="velocity_disp_theta")
        svnm_theta_E   = save_json_name(setting=sett,path=svpth,filename="theta_E_mcmc")
        svnm_M_theta_E = save_json_name(setting=sett,path=svpth,filename="M_theta")
        svnm_M_E_mcmc  = save_json_name(setting=sett,path=svpth,filename="M_Einst")
        svnm_m_cumul   = save_json_name(setting=sett,path=svpth,filename="m_cumul")
        if prior_res:  
            try:
                sig_theta_mcmc = load_whatever(svnm_sig_theta)*u.km/u.second
                theta_E_mcmc   = load_whatever(svnm_theta_E)*u.arcsec
                M_theta_mcmc   = load_whatever(svnm_M_theta_E)*u.solMass
                M_E_mcmc       = load_whatever(svnm_M_E_mcmc)*u.solMass
                m_cumul        = load_whatever(svnm_m_cumul)*u.solMass
                print("Results loaded")
            except:
                print("No previous results found")
                prior_res = False
        if not prior_res:
            print("Calculating...")
            sig_theta_mcmc,M_E_mcmc,M_theta_mcmc,m_cumul,theta_E_mcmc = get_M_w_R(sett,ret_w_units=True,cut_mcmc=cut_mcmc)
            
            save_json(sig_theta_mcmc.value,svnm_sig_theta)
            save_json(theta_E_mcmc.value,svnm_theta_E)
            save_json(M_theta_mcmc.value,svnm_M_theta_E)
            save_json(M_E_mcmc.value,svnm_M_E_mcmc)
            save_json(m_cumul.value,svnm_m_cumul)
        sigtheta_res_str   = easy_res_str(sig_theta_mcmc.value,quantiles=True)
        print("sigma_theta =",sigtheta_res_str,sig_theta_mcmc.unit)
        M_E_mcmc_res_str   = easy_res_str(M_E_mcmc.value,quantiles=True)
        print("M_Einstein  =",M_E_mcmc_res_str,M_E_mcmc.unit)
        M_theta_mc_res_str = easy_res_str(M_theta_mcmc.value,quantiles=True)
        print("M(theta) =",M_theta_mc_res_str,M_theta_mcmc.unit)

        #######
        # MtL #
        #######
        #m_cumul shape = n_mcmc,r_i
        mass_cumul_mcmc = np.transpose(m_cumul)
        mass_sig_up,mass,mass_sig_low = [],[],[]
        for i,ri in enumerate(sett.rad):
            #mass_cumul_mcmc shape: r_i, n_mcmc
            min_m,med_m,max_m = quantile(mass_cumul_mcmc[i],q=[0.16,0.5,0.84])
            try:
                mass.append(med_m.value)
                mass_sig_low.append(min_m.value)
                mass_sig_up.append(max_m.value)            
            except AttributeError:
                mass.append(med_m)
                mass_sig_low.append(min_m)
                mass_sig_up.append(max_m)
        mass         = np.array(mass)*u.solMass
        mass_sig_low = np.array(mass_sig_low)*u.solMass
        mass_sig_up  = np.array(mass_sig_up)*u.solMass
        
        mtl = np.array(mass)/np.array(sett.Lum)
        index_cutcenter = np.where(mtl==np.min(mtl))[0][0]

        fig, ax = plt.subplots(3,figsize=(9,12))
        
        print(r"$r_{min}$"+f" ={sett.rad[index_cutcenter]}")
        print(r"$r_{min}$"+f" ={sett.rad[index_cutcenter]/cosmo.arcsec_per_kpc_proper(sett.z_lens)}")
        print(r"$r_{min}[pix]$"+f" ={sett.rad[index_cutcenter]/sett.pix_scale}")
        theta_E_mcmc                 = theta_E_mcmc.value
        min_thetaE,thetaE,max_thetaE = quantile(theta_E_mcmc,q=[.16,0.5,.84])
        res_thetaE      = easy_res_str(theta_E_mcmc,quantiles=True) 
        index_theta_E   = int(np.where(np.abs(np.median(theta_E_mcmc)-sett.rad.value[index_cutcenter:])==np.min(np.abs(np.median(theta_E_mcmc)-sett.rad.value[index_cutcenter:])))[0])
        N_thetaE_cut    = 10
        tomin           = np.abs((np.median(theta_E_mcmc)*N_thetaE_cut)-sett.rad.value[index_cutcenter:])
        index_N_theta_E = int(np.where(tomin==np.min(tomin))[0])+index_cutcenter
        cut_theta_E     = True
        if not cut_theta_E:
            index_N_theta_E = len(mtl)
        else:
            print("Cutting plots at ",N_thetaE_cut," time the theta_E, equals to ",np.round(sett.rad[index_N_theta_E],2))
        mtl  = mtl[index_cutcenter:index_N_theta_E]
        mass = mass[index_cutcenter:index_N_theta_E]
        lum  = np.array(sett.Lum)[index_cutcenter:index_N_theta_E]
        rad  = sett.rad[index_cutcenter:index_N_theta_E]
        min_mass = np.array(mass_sig_low)[index_cutcenter:index_N_theta_E]
        max_mass = np.array(mass_sig_up)[index_cutcenter:index_N_theta_E]
        min_mtl  = min_mass/lum
        max_mtl  = max_mass/lum 
        ax[0].set_title(r"Mass [$M_\odot$]")
        ax[0].fill_between(rad,max_mass,min_mass,color="cyan",alpha=.5,label=r"M 1-$\sigma$ region")
        ax[0].plot(rad,mass,color="k")
        mass = mass.value
        res_mass = simp(mass[index_theta_E],max_mass[index_theta_E]-mass[index_theta_E],mass[index_theta_E]-min_mass[index_theta_E])+r"$M_\odot$"
        ax[0].plot(rad,mass,color="k",alpha=0,label=r"M(<$\theta_E$>)= "+res_mass)
        ax[0].set_ylabel(r'M(a) [M$_\odot$]') 
        
        ax[1].set_title(r"Luminosity [$L_\odot$]")
        ax[1].plot(rad,lum,color="k")
        err_lum = np.array(sett.errLum)[index_cutcenter:index_N_theta_E]
        max_lum = lum+err_lum
        min_lum = lum-err_lum
        ax[1].fill_between(rad,max_lum,min_lum,color="cyan",alpha=.5,label=r"L 1-$\sigma$ region")
        res_lum = simp(lum[index_theta_E],err_lum[index_theta_E],err_lum[index_theta_E])
        ax[1].plot(rad,lum,color="k",alpha=0,label=r"L(<$\theta_E$>)= "+res_lum+r"$L_\odot$")
        ax[1].set_ylabel(r'L(a) [$L_\odot$]') 

        ax[2].set_title(r"Mass to Light ratio $\mathcal{\Upsilon}$")
        ax[2].fill_between(rad,max_mtl,min_mtl,color="cyan",alpha=.5,label=r"$\mathcal{\Upsilon}$ 1-$\sigma$ region")
        ax[2].plot(rad,mtl,color="k")
        res_mtl = simp(mtl[index_theta_E],max_mtl[index_theta_E]-mtl[index_theta_E],mtl[index_theta_E]-min_mtl[index_theta_E])+r"$M_\odot$/$L_\odot$"
        ax[2].plot(rad,lum,color="k",alpha=0,label=r"$\mathcal{\Upsilon}$(<$\theta_E$>)= "+res_mtl)
        ax[2].set_ylabel(r'$\mathcal{\Upsilon}$(a) [$M_\odot$/$L_\odot$]')

        for i,axi in enumerate(ax): 
            if i==1:
                axi.axvspan(min_thetaE,max_thetaE, color="grey",alpha=.5,label=r"$\theta_E$ 1-$\sigma$ region")
                axi.axvline(np.median(theta_E_mcmc),c="r",ls="--",label=r"<$\theta_E$>="+res_thetaE+" ''")
            else:
                axi.axvspan(min_thetaE,max_thetaE, color="grey",alpha=.5)#,label=r"$\theta_E$ 1-$\sigma$ region")
                axi.axvline(np.median(theta_E_mcmc),c="r",ls="--")#,label=r"<$\theta_E$>="+res_thetaE+" ''")
            axi.set_xlabel('a [arcsec]') #it's not r, it's the semi-major axis
            ax_tmp = axi.twiny() 
            ax_tmp.set_xlabel('a [kpc]')
            ax_tmp.plot(rad/cosmo.arcsec_per_kpc_proper(sett.z_lens), np.zeros_like(rad),alpha=0.)
            ax_tmp.tick_params(axis='x')
            axi.legend(loc='lower right')
            if i!=len(ax)-1:
                oft = axi.yaxis.get_offset_text()
                oft.set_x(-0.04)
        # for some reason the axis ratio is the same for all, which is not good for the mtl
        maxmaxmtl,minminmtl = np.max(max_mtl),np.min(min_mtl)
        dmtl = (maxmaxmtl-minminmtl)*1./10
        ax[2].set_ylim(minminmtl-dmtl,maxmaxmtl+dmtl)
        mtl_name = f"{get_savefigpath(sett)}/MtL_rewamp_2.pdf"
        print(f"Saved {mtl_name}")
        plt.tight_layout()
        plt.savefig(mtl_name)
        plt.close()

    produce_color_plots(sett1=sett_f475,sett2=sett_f814,max_rad=3)
    
    success(sys.argv[0])
