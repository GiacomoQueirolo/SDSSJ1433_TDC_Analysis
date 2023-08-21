#models = ["isopy_f475w_220525.fits","model_tot_f814w_220527.fits"] #,"new_outdata_model_f105w.fits","new2_outdata_model_f140w.fits","model_tot_f160w_030423.fits"]

import sys
import argparse
import numpy as np
from acstools import acszpt 
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy import units as unit
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

from Utils.tools import *
from Utils.get_res import *
from Data.conversion import *

from Data.input_data import init_lens_model
from Data.image_manipulation import load_fits,get_header

@check_setting
def getluminosity(setting,intens,r,cosmo=None,K=None):
    if cosmo is None:
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # cosmo from https://academic.oup.com/mnras/article/474/3/3391/4644836, Agnello 2017
    if K is None:
            K = get_K_corr(setting=setting)
    cosmo_ld = cosmo.luminosity_distance(z=setting.z_lens).to("pc")
    filter   = get_filter(setting)
    ZP       = get_ZP(filt=filter)
    M_sun    = get_M_sun(filt=filter)
    return _getluminosity(intens,M_sun=M_sun,r=r,distance=cosmo_ld,ZP=ZP,pix_scale=setting.pix_scale,K=K)

def get_K_corr(setting):
    filter = get_filter(setting)
    # hardcoded from http://kcor.sai.msu.ru/
    # z_lens = 0.407
    # Johnson B, color B-Ic =1.2
    if "475" in filter:
        return 0.09 
    # Cousins Ic , color V-Ic =0.5
    elif "814" in filter:
        return -0.14
    else:
        #pragma no cover
        raise RuntimeError("Not implemented")

 
def _getluminosity(intens,M_sun,r,eps,distance,ZP,pix_scale,z_lens,K):

    flux=[]
    area_annulus = get_annulus(r=r,eps=eps,from_0=True)
    flux_i = area_annulus*intens[i]
    flux =[flux_i[0]]
    for i in range(len(intens)):
        if i==0:
            flux.append(flux_i[0])
        else:
            flux.append(flux_i[i]+flux[-1])
    # x is in arcsec
    for i in range(len(intens)):
        if i==0:
            flux.append(np.pi*(r[i]**2)*intens[i])
        else:
            flux.append(flux[i-1] + np.pi*(r[i]**2-r[i-1]**2)*intens[i])

    mag =[] # AB mag
    Mag =[] # AB mag
    Lum =[] # solar_lum

    d_mod=5*np.log10(distance.to("pc")/(10*unit.pc))
    cosmic_dimming = 2.5*np.log10((1+z_lens)**4)

    for i in range(len(flux)):
        mag.append(-2.5*np.log10(flux[i]/(pix_scale**2))+ZP)
        Mag.append(mag[i]-d_mod-K -cosmic_dimming )#consider galactic absorption , cosmic dimming
        Lum.append(10**(0.4*(M_sun-Mag[i])))

    return Lum,Mag,mag

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

def velocity_disp_FJr(L_V):
    # unit L_v = solar Lum (total)
    return (200*unit.km/unit.s)*((L_V/(2*10**10))**(1/4))


def velocity_disp_theta(theta,z_lens,z_source,cosmo=None):
    if cosmo is None:
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    cosmo_ds  = cosmo.angular_diameter_distance(z_source)
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source)
    theta_rad = np.pi*theta/(180*3600)
    return const.c.to("km/s")*np.sqrt(theta_rad*cosmo_ds/(4*np.pi*cosmo_dds))

def get_M_sun(filt):
    if "475" in filt:
        M_sun=5.07# AB_mag
    elif "814" in filt:
        M_sun=4.52 #AB_mag 
    else:
        raise RuntimeError(f"No M_sun defined for filter {filt}")
    return M_sun
    
def spline_lum(x_fit,x_value,y_value):
    tck = interpolate.splrep(x_value, y_value)
    return interpolate.splev(x_fit, tck)

def get_kappa_along_r(r,kwres,lens_model):
    # have to be centered on the lens center
    x_0   = kwres["kwargs_lens"][0]["center_x"]
    y_0   = kwres["kwargs_lens"][0]["center_y"]    
    x0    = x_0*np.ones_like(r)
    y_r   = y_0+np.array(r)
    # Test
    ###########################
    """
    # Seems fine to me
    if sett is not None:
        from reg_mask import str_reg
        #x_0_px,y_0_px = conv_radec_to_xy(sett,ra=x_0,dec=y_0)
        #for ri in r:
            #str_reg+=f"circle({x_0_px[0]+1},{y_0_px[0]+1},{ri/sett.pix_scale}) # color=chocolate dash=1\n"
        x_px,y_px = conv_radec_to_xy(sett,ra=x0,dec=y_r)
        for xi,yi in zip(x_px,y_px):
            str_reg+=f"circle({xi+1},{yi+1},3) # color=green dash=1\n"
        reg_file = get_savefigpath(setting=sett)+"/_tmp_kappa_rad.reg"
        with open(reg_file, "w") as text_file:
                print(str_reg, file=text_file)
        print("TEST: saved ",reg_file)
    """
    ###########################
    kappa = lens_model.kappa(x=x0,y=y_r,kwargs=kwres['kwargs_lens'])
    return kappa

from Data.image_manipulation import get_rotangle
def get_kappa(r,eps,pa,kwres,lens_model,transform_pix2angle,phi_thresh=5,q_thresh=.1):
    # have to be centered on the lens center
    x_0   = kwres["kwargs_lens"][0]["center_x"]
    y_0   = kwres["kwargs_lens"][0]["center_y"]    
    # 2 options: A) ~ where we ignore the difference of pointing btw lens and light
    # B) we consider the difference and integrate over the annulus
    # for now A)
    
    # r is assumed to be in arcsecs
    # pa is assumed to be in deg 
    rotang = get_rotangle(transform_pix2angle)
    phi    = (rotang-np.array(pa)) - 180
    if np.mean(phi)<-100:
        phi+=180
    phi_rad = np.pi*phi/(180*3600)
    x = x_0 + np.array(r)*np.cos(phi_rad) #x_0*np.ones_like(r)
    y = y_0 + np.array(r)*np.sin(phi_rad) #x_0*np.ones_like(r)
    
    # B) for now only test/warning
    q  = 1 - np.array(eps)
    e1,e2 = kwres["kwargs_lens"][0]["e1"],kwres["kwargs_lens"][0]["e2"]
    q_lm,phi_lm = qphi_from_e1e2(e1,e2)
    if np.abs(phi-phi_lm)>phi_thresh:
        print(f"WARNING: Difference of angle higher then {phi_thresh}: {np.tot(np.abs(phi-phi_lm))},")
    if np.abs(q-q_lm)>q_thresh:
        print(f"WARNING: Difference of ellipticity higher then {q_thresh}: {np.tot(np.abs(q-q_lm))},")
    kappa = lens_model.kappa(x=x,y=y,kwargs=kwres['kwargs_lens'])
    return kappa

def get_annulus(r,eps,from_0=True):
    # r is considered in arcsec or px, not **1/4
    # copied from sersic_def.py of Matthias Kluge
    area_ann = np.pi*(r[1:]**2*(1-eps[1:]) - r[:-1]**2*(1-eps[:-1]) )
    if from_0:
        list_ann = np.array(area_ann).tolist()
        area_ann = [np.pi*r[0]**2*(1-eps[0]), *list_ann ]
    return area_ann

def _get_cumulative_mass(kappa,r,crit_density):
    cumulative_mass=[]
    if getattr(crit_density,"value",False):
        crit_density = crit_density.to("Msun/(arcsec arcsec)").value
        print("Crit density ",crit_density/1e9,"1e9 Msun/arcsec^2")
    cumulative_mass = [np.pi*(r[0])**2*kappa[0]*crit_density]
    for i in range(1,len(r)):
        cumulative_mass.append(np.pi*(r[i]**2-r[i-1]**2)*kappa[i]*crit_density+cumulative_mass[i-1])
    return cumulative_mass

@check_setting
def get_cosmo_prm(setting,cosmo=None,SigCr_arcs2=False):
    # SigCr_arcs2: if True, convert Sigma_Crit in M_sun/arcsec^2
    # instead of M_sun/kpc^2
    z_source  = setting.z_source
    z_lens    = setting.z_lens
    if not cosmo:
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # cosmo from https://academic.oup.com/mnras/article/474/3/3391/4644836, Agnello 2017
    cosmo_dd  = cosmo.angular_diameter_distance(z_lens).to("kpc")   #kpc
    cosmo_ds  = cosmo.angular_diameter_distance(z_source).to("kpc") #kpc
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source).to("kpc") #kpc
    
    # Sigma_crit = D_s*c^2/(4piG * D_d*D_ds)
    Sigma_crit   = (cosmo_ds*const.c**2)/(4*np.pi*const.G*cosmo_dds*cosmo_dd)
    Sigma_crit   = Sigma_crit.to("Msun /(kpc kpc)")    
    if SigCr_arcs2:
        Sigma_crit/=(cosmo.arcsec_per_kpc_proper(z_lens))**2
    return Sigma_crit,cosmo_dd

def get_cumulative_mass(sett,r,eps,pa,kwres=None,lens_model=None,cosmo=None):
    if lens_model is None:
        lens_model = init_lens_model(sett)
    if kwres is None:
        kwres = get_kwres(sett)["kwargs_results"]
    if cosmo is None:
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # cosmo from https://academic.oup.com/mnras/article/474/3/3391/4644836, Agnello 2017
    #kappa           = get_kappa_along_r(r,kwres=kwres,lens_model=lens_model)
    kappa          = get_kappa(r,eps,pa,kwres=kwres,lens_model=lens_model,transform_pix2angle=sett.transform_pix2angle)
    Sigma_crit, _,_ = get_cosmo_prm(sett,cosmo=cosmo,SigCr_arcs2=True)
    cumulative_mass = _get_cumulative_mass(kappa=kappa,r=r,crit_density=Sigma_crit)
    return cumulative_mass

def get_radii_images(setting,kw_res=None):
    # Consider the mass enclosed in a radius of approximately at the image positions, 
    # centered on the center of the main lens
    if kw_res is None:
        kw_res = get_kwres(setting)
    lens_res            = kw_res["kwargs_lens"][0]
    x,y                 = lens_res["center_x"],lens_res["center_y"]
    ps_res              = kw_res["kwargs_ps"][0]
    ra_image, dec_image = ps_res["ra_image"],ps_res["dec_image"] # A is ~0,0
    radii  = [np.sqrt((ra_i-x)**2 + (dec_i-y)**2) for ra_i,dec_i in zip(ra_image,dec_image)]
    return radii 

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(prog=sys.argv[0],description="Mass to light ratio",formatter_class=CustomFormatter)
    #parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args  = parser.parse_args()
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # cosmo from https://academic.oup.com/mnras/article/474/3/3391/4644836, Agnello 2017
    sett_f475 = get_setting_module("f475w_CPRP_PLLI_ws",1)
    sett_f814 = get_setting_module("f814w_CPRP_PLLI_ws",1)
    for sett in [sett_f475,sett_f814]:
        data_path  = sett.data_path
        model_name = sett.lens_light_model_name_full
        if "814" in get_filter(sett):
            model_name = "isopy_220527.fits"
        model_path = f"{data_path}/{model_name}"

        param_val = load_fits(model_path,HDU=-1)
        a14       = param_val.field("a14")[1:]
        sett.rad  = a14**4 #arcsec
        intens    = param_val.field("flux_sb")[1:]
        sett.eps_light = param_val.field("eps")[1:]
        sett.pa_light  = param_val.field("pa")[1:]



        sett.L,sett.Mag,sett.mag  = getluminosity(sett,intens,r=sett.rad,cosmo=cosmo,K=0)
        ########################################################################
        # SB is also calculated, but the ZP is wrong
        # correct for that and check that it's compatible with what we obtain
        ########################################################################
        SB_tc    = param_val.field("sb")[1:]
        try:
            ZP_head = get_header(f"{data_path}/{sett.image_name}","ZP")
        except KeyError:
            print("ZP in header not found, assumed to be 0")
            ZP_head = 0
        ZP_litt  = get_ZP(get_filter(sett))
        DZP      = ZP_head - ZP_litt
        SB       = np.array(SB_tc)-DZP
        area     = np.pi*np.array(sett.rad)**2
        _mag     = ZP_litt -2.5*np.log10(np.array(intens)*area/(sett.pix_scale**2))
        SB_from_mag = _mag +2.5*np.log10(area)
        print("Test SB")
        D_SB = (SB-SB_from_mag)
        mean_dsb,std_dsb = np.nanmean(D_SB),np.nanstd(D_SB)
        if abs(mean_dsb)>1e-3 or std_dsb>1e-10: 
            print("Surface brightness not corrected")
            print("mean",mean_dsb)
            print("std",std_dsb)
            D_SB = D_SB[~np.isnan(D_SB)]
            plt.close()
            plt.hist(np.array(D_SB).tolist())
            plt.savefig(get_savefigpath(sett)+"/.D_SB.pdf")
            print("Saved "+get_savefigpath(sett)+"/.D_SB.pdf")
            plt.close()
        else:
            print("Test passed")
        #SB_cumul =  sett.mag +2.5*np.log10(area)
        plt.plot(sett.rad,sett.Mag )
        plt.xlabel("rad [arcsec]")
        plt.ylabel(r"Abs. Mag [mag]")
        plt.savefig(get_savefigpath(sett)+"/Mag.pdf")
        print("Saved "+get_savefigpath(sett)+"/Mag.pdf")
        plt.close() 
        ########################################################################

    """
    sett_f814.mag_interp = spline_lum(x_fit=sett_f475.rad,x_value=sett_f814.rad,y_value=sett_f814.mag)
    color_f475_f814 = np.array(sett_f475.mag)-np.array(sett_f814.mag_interp)
    plt.plot(sett_f475.rad,color_f475_f814,'.')
    plt.xlabel('r [arcsec]')
    plt.ylabel('Apparent F475-F814')
    #plt.xlim([0,1.4])
    for sett in [sett_f475,sett_f814]:
        color_plot_name = f"{get_savefigpath(sett)}/apparent_color_f475_f814.pdf"
        plt.savefig(color_plot_name)
        print(f"Saved {color_plot_name}") 
    plt.close()
    new_rad = np.linspace(0,9,num=100)
    mag_interp_f814 = spline_lum(x_fit=new_rad,x_value=sett_f814.rad,y_value=sett_f814.mag)
    mag_interp_f475 = spline_lum(x_fit=new_rad,x_value=sett_f475.rad,y_value=sett_f475.mag)
    color_f475_f814_interp = np.array(mag_interp_f475)-np.array(mag_interp_f814)
    plt.plot(new_rad,color_f475_f814_interp,'.')
    plt.xlabel('r [arcsec]')
    plt.ylabel('Apparent F475-F814 interpolated')
    #plt.xlim([0,1.4])
    for sett in [sett_f475,sett_f814]:
        color_plot_name = f"{get_savefigpath(sett)}/apparent_color_f475_f814_interpolated.pdf"
        plt.savefig(color_plot_name)
        print(f"Saved {color_plot_name}") 
    plt.close()
    """
    new_rad = np.linspace(0,9,num=100)
    mag_interp_f814 = spline_lum(x_fit=new_rad,x_value=sett_f814.rad,y_value=sett_f814.mag)
    mag_interp_f475 = spline_lum(x_fit=new_rad,x_value=sett_f475.rad,y_value=sett_f475.mag)

    Mag_interp_f814 = spline_lum(x_fit=new_rad,x_value=sett_f814.rad,y_value=sett_f814.Mag)
    Mag_interp_f475 = spline_lum(x_fit=new_rad,x_value=sett_f475.rad,y_value=sett_f475.Mag)
    abs_color_f475_f814_interp = np.array(Mag_interp_f475)-np.array(Mag_interp_f814)
    plt.plot(new_rad,abs_color_f475_f814_interp,'.')
    plt.xlabel('r [arcsec]')
    plt.ylabel('Absolute F475-F814 interpolated')
    plt.title("Absolute color")
    for sett in [sett_f475,sett_f814]:
        color_plot_name = f"{get_savefigpath(sett)}/absolute_color_f475_f814_interpolated.pdf"
        plt.savefig(color_plot_name)
        print(f"Saved {color_plot_name}") 
    plt.close()
    
    plt.plot(new_rad, mag_interp_f814,".",label="f814w")
    plt.plot(new_rad, mag_interp_f475,".",label="f475w")
    plt.legend()
    plt.xlabel('r [arcsec]')
    plt.ylabel('mag') 
    for sett in [sett_f475,sett_f814]:
        color_plot_name = f"{get_savefigpath(sett)}/mag_both_filters.pdf"
        plt.savefig(color_plot_name)
        print(f"Saved {color_plot_name}") 
    plt.close()

    ######################################
    # Mass and MtL ratio
    lens_model = init_lens_model(setting=sett_f814) #should be the same for both sett
    for sett in [sett_f475,sett_f814]:
        sett.kwres      = get_kwres(sett)["kwargs_results"]
        print("sigma_FJr =",velocity_disp_FJr(sett.Lum[-1]))
        sigma_theta = velocity_disp_theta(sett.kwres["kwargs_lens"][0]["theta_E"],sett.z_lens,sett.z_source,cosmo)
        print("sigma_theta =",sigma_theta)

        sett.mass_cumul = get_cumulative_mass(sett=sett,r=sett.rad,eps=sett.eps_light,pa=sett.pa_light,kwres=sett.kwres,\
                                              lens_model=lens_model,cosmo=cosmo)
        fig, ax1 = plt.subplots()
        ax1.set_title("Cumulative Mass")
        ax1.plot(sett.rad, sett.mass_cumul,".",alpha=0)
        ax1.set_xlabel('r [arcsec]')
        ax1.set_ylabel(r'Mass cumulative [M$_\odot$]')
        ax2 = ax1.twiny() 
        ax2.set_xlabel('r [kpc]')
        ax2.plot(sett.rad/cosmo.arcsec_per_kpc_proper(sett.z_lens), sett.mass_cumul,".")
        ax2.tick_params(axis='x')
        plot_name = f"{get_savefigpath(sett)}/Mass_cumul.pdf"
        plt.tight_layout()
        plt.savefig(plot_name)
        print(f"Saved {plot_name}") 
        plt.close()

        fig, ax1 = plt.subplots()
        ax1.set_title("Cumulative Luminosity")
        ax1.plot(sett.rad, sett.Lum,".",alpha=0)
        ax1.set_xlabel('r [arcsec]')
        ax1.set_ylabel(r'Lum cumulative [Lum$_\odot$]')      
        ax2 = ax1.twiny() 
        ax2.set_xlabel('r [kpc]')
        ax2.plot(sett.rad/cosmo.arcsec_per_kpc_proper(sett.z_lens), sett.Lum,".")
        ax2.tick_params(axis='x')
        plot_name = f"{get_savefigpath(sett)}/L_cumul.pdf"
        plt.tight_layout()
        plt.savefig(plot_name)
        print(f"Saved {plot_name}") 
        plt.close()
        
        mtl = np.array(sett.mass_cumul)/np.array(sett.Lum) # M_sun/L_sun
        fig, ax1 = plt.subplots()
        ax1.set_title("Mass-to-Light Ratio")
        ax1.plot(sett.rad, mtl,alpha=0,label=get_filter(sett))
        ax1.axvline(get_radii_images(sett,sett.kwres)[0],c="k",label="Image A")
        ax1.set_xlabel('r [arcsec]')
        ax1.set_ylabel(r'M/Lum [M$_\odot$/Lum$_\odot$]')
        ax1.legend()
        ax2 = ax1.twiny() 
        ax2.set_xlabel('r [kpc]')
        ax2.plot(sett.rad/cosmo.arcsec_per_kpc_proper(sett.z_lens), mtl,label=get_filter(sett))
        ax2.tick_params(axis='x')
        mtl_name = f"{get_savefigpath(sett)}/MtL_rewamp.pdf"
        plt.tight_layout()
        plt.savefig(mtl_name)
        print(f"Saved {mtl_name}") 
        plt.close()        
    success(sys.argv[0])