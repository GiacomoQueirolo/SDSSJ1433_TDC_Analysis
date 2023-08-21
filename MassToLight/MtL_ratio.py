from Utils.tools import *
from Utils.get_res import *
from Data.image_manipulation import *
from MassToLight.grid_class import Circ_Grid
from MassToLight.pos_img import create_pos_img
from Utils.create_fits import lens_light
from Data.input_data import init_lens_model
from Data.conversion import conv_radec_to_xy,conv_xy_to_radec

import sys
import pickle
import argparse
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy.cosmology import default_cosmology
   
def get_cosmo_prm(setting,cosmo=None):
    z_source  = setting.z_source
    z_lens    = setting.z_lens
    if not cosmo:
        cosmo = default_cosmology.get()
    cosmo_dd  = cosmo.angular_diameter_distance(z_lens).to("kpc")   #kpc
    cosmo_ds  = cosmo.angular_diameter_distance(z_source).to("kpc") #kpc
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source).to("kpc") #kpc
    
    # Sigma_crit = c^2*D_s/(4piG * D_d*D_ds)
    Sigma_crit = const.c*const.c*cosmo_ds/(4*np.pi*const.G*cosmo_dds*cosmo_dd)
    Sigma_crit = Sigma_crit.to("Msun/(kpc kpc)") # = to the sigmacrit of lenstronomy up to 0.03%
    
    cosmo_Ld  = cosmo.luminosity_distance(z_lens) 
    return Sigma_crit,cosmo_dd,cosmo_Ld

def arcsec_to_kpc(dx,cosmoDd,rounded=False):
    kpc = dx*u.arcsec.to("rad")*cosmoDd.to("kpc")
    if rounded:
        return format(kpc,".3f")
    return kpc


def get_kappa_discreet(setting,radius,lens_model,kwlens,rmask,dr=0.1,input_type="radec",grid=None):
    center_lens = kwlens[0]["center_x"],kwlens[0]["center_y"] 
    if grid is None:
        grid    = Circ_Grid(center=center_lens,radius=radius,edge=dr,input_type=input_type,setting=setting,r_mask=rmask)
    grid_circ   = grid.get_circ_grid()
    gcf         = np.transpose([gcf.get_radec() for gcf in grid_circ]).tolist() 
    kappa_flat  = lens_model.kappa(*gcf,kwargs=kwlens)
    return kappa_flat,grid


def get_mass_grid(setting,radius,lens_model,kwlens,Sigma_crit,cosmo_dd,rmask,dr=0.1,input_type="radec",grid=None):
    # dr= precision in arcsec
    # Mass(r) = 2 pi Sigma_crit * integr_0^r dr kappa(r) * r
    # or rather, since it's discretised anyway:
    # density(i) = Sigma_cir * kappa(i) with i indicating the grid cells within the radius
    # area(i)[Mpc] = area(i)[rad]*Dd**2
    # M(i) = density(i)*area(i)[kpc]
    # M(r) = Sum_[i<=r] M(i)
    # equivalent to:  
    # M(r) = Sigma_crit * Sum_i kappa(i)*Area(i) 
    kappa_discreet, grid = get_kappa_discreet(setting,radius,lens_model,kwlens,rmask,dr=dr,input_type=input_type,grid=grid) 
    density_grid         = kappa_discreet*Sigma_crit # M_Sun/kpc^2
    area_pix             = grid.get_area_square()*np.ones_like(grid.get_circ_grid())
    #area_pix_kpc2        = area_pix*(u.arcsec.to("rad")*cosmo.arcsec_per_kpc_proper(z_lens)**2)#cosmo_dd*cosmo_dd*area_pix*u.arcsec.to("rad")*u.arcsec.to("rad") #
    area_pix_kpc2        = cosmo_dd*cosmo_dd*area_pix*u.arcsec.to("rad")*u.arcsec.to("rad") 
    mass_grid            = density_grid*area_pix_kpc2  # M_Sun
    #Mass = np.sum(mass_grid)
    #Mass = Mass.to("Msun")
    #mass_3D   = lensModel.lens_model.mass_3d(r=radius,kwargs=kwlens["kwargs_lens"],bool_list=[0]) #bc shear do not have the mass and SIS is ignored (see note may 23th '22)
    #mass_fact = u.arcsec.to("rad") ** 2 * cosmo_dd * cosmo_ds / cosmo_dds * const.c ** 2 / (4 * np.pi * const.G)

    #mass = mass_3D*mass_fact
    #mass = mass.to("Msun") #solar masses
    #return Mass,mass_grid,grid 
    return mass_grid,grid 

def plot_MtR(setting,mass_grid,grid,savename=None):
    dists = [gp.dist_from_center.length_radec for gp in grid.get_circ_grid()]

    r0,r1 = min(dists),max(dists)
    dr    = grid.edge_radec
    bins  = np.arange(r0,r1,dr)
    sorted_indexes = np.argsort(dists).tolist()
    sorted_mass    = np.array(mass_grid)[sorted_indexes]
    sorted_dist    = np.array(dists)[sorted_indexes]
    m = np.zeros(shape=len(bins))
    r = bins+dr/2.
    bin0 = np.where(sorted_dist<=min(bins))
    m[0] = np.sum(sorted_mass[bin0])
    sorted_dist = sorted_dist[bin0[0][-1]:]
    sorted_mass = sorted_mass[bin0[0][-1]:]
    
    for i,d in enumerate(sorted_dist):
        if d<=min(bins):
            m[0]+=sorted_mass[i]
        else:
            index_bin = np.array(np.where((bins>d) & (bins<d+dr)))-1
            for ii in index_bin:
                m[ii] += sorted_mass[i]
    plt.scatter(r,m,marker="x",color="k")
    plt.xlabel("Distance from Lens Center [\"]")
    plt.ylabel("Mass [Msun]")
    plt.title("M(R)")
    if savename is not None:
        plt.savefig(savename)
        print("Saved plot "+savename)
    plt.close()
    return 0 

#################################
# get light from pixellated map #
#################################

def flux_in_grid(lens_light_pix,grid): 
    # now the problem is, the grid is not discretised in pixels
    grid_degraded = grid.get_degraded_pixel_grid()
    # we now have a grid with the coord of the pixels
    flux = np.sum([lens_light_pix[i[0]][i[1]] for i in grid_degraded])
    return flux
    
if __name__=="__main__":
    present_program(sys.argv[0])

    parser = argparse.ArgumentParser(description="Mass, Light and Mass-to-Light ratio within a certain radius for given models")
    #parser.add_argument("-PR", "--plot_wrt_radius",dest="plot_wrt_radius", default=False,action="store_true",
    #                    help="Plot with respect to radius")
    parser.add_argument("-r", "--radius", type=float, dest="radius", default=-1,
                        help="Radius in arcsec at which to calculate the mass. Default: approx max radius of the lensed image position")                 
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting files to consider")
    
    args      = parser.parse_args()
    #PR        = args.plot_wrt_radius
    radius    = args.radius
    setting_names  = args.SETTING_FILES
    ###
    dr    = 0.005 #arcsec
    rmask = 0.3   #arcsec
    ###
    to_redo = False

    for sett in setting_names:
        print_setting(sett)
        setting   = get_setting_module(sett,1)
        kw_res    = get_kwres(setting)["kwargs_results"]
        # we ignore the perturber and the shear
        kwlens    = [kw_res["kwargs_lens"][0]]
        lensModel = init_lens_model(setting,["PEMD"])
        SigmaCrit,cosmoDd,cosmoLd = get_cosmo_prm(setting)
        
        print("dr :",arcsec_to_kpc(dr,cosmoDd,rounded=True))
        if radius==-1 or to_redo:
            print("Consider the mass enclosed in a radius of approximately at the image positions, centered on the center of the main lens")
            #print("from lenstronomy: \"This function assumes spherical symmetry/ignoring the eccentricity.\" ")
            lens_res            = kw_res["kwargs_lens"][0]
            x,y                 = lens_res["center_x"],lens_res["center_y"]
            ps_res              = kw_res["kwargs_ps"][0]
            ra_image, dec_image = ps_res["ra_image"],ps_res["dec_image"] # A is ~0,0
            max_darcsec = [max([np.sqrt((ra_i-x)**2 + (dec_i-y)**2) for ra_i,dec_i in zip(ra_image,dec_image)])]
            radius  = max(max_darcsec)/2.
            to_redo = True
        grid_path = get_savefigpath(setting)+"/grid_MtL_dr"+str(dr)+"_r"+str(radius)+".pkl"
        try:
            grid = pickle.load(open(grid_path,"rb"))
        except:
            grid = None 
        mass_grid,grid = get_mass_grid(setting,radius,lensModel,kwlens,SigmaCrit,cosmoDd,rmask=rmask,dr=dr,input_type="radec",grid=grid)
        mass_discreet = np.sum(mass_grid).to("Msun") # solar masses
        plot_MtR(setting,mass_grid,grid,savename=get_savefigpath(setting)+"/MtR.pdf")
        print("Radius ", format(radius,".2f"),"\"= ",arcsec_to_kpc(radius,cosmoDd,True))
        print("Mass enclosed in radius ", format(radius,".2f"),"\": ",format(mass_discreet.value,".2E")," Msun")
        # plot pixellated grid: 
        # TO TEST
        with open(grid_path,"wb") as f:
            pickle.dump(grid,f)
        print("Saved grid as pickle in "+grid_path)
        ax_im = create_pos_img(setting,RP=False,RM=False,save=False,also_fits=False)
        ax_im_pix = grid.plot_grid(ax=ax_im, grid_to_plot="circ_grid")
        plt.savefig(get_savefigpath(setting)+"/grid_plot.pdf")
        plt.close()
        grid.produce_degraded_fits()
        flux_unit    = 1/u.second #count per sec
        #https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints
        #https://hst-docs.stsci.edu/wfc3dhb/chapter-9-wfc3-data-analysis/9-1-photometry
        if check_if_SUB(sett): # this is SUB=> model from isophotes
            lens_light_pix = get_lens_light_model(sett)
            path_llm = setting.data_path+"/"+setting.lens_light_model_name
            #ZP_mag   = get_header(path_llm,"ZP")
            #photonu  = get_header(path_llm,"PHOTFNU")*u.second/1 # ~Jy~ * sec/count -> not Jy ? must somehow already be taken into account  
            #if "475" in get_filter(sett) or "814" in get_filter(sett):
            #   UVIS : https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/uvis-encircled-energy
            #   IR : https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/ir-photometric-calibration
            #   for both Flux = Fr*photoflam/EE(r)
            wvlength  = get_header(path_llm,"PHOTPLAM")# Angstrom
            bandwidth = get_header(path_llm,"PHOTBW")# Angstrom
            PHOTFLAM  = get_PHOTFLAM(setting)
            EE_r      = get_EE(setting,radius)
            flux_fact = PHOTFLAM/EE_r # ergs/cm2/Angs/e-
            # should this be multiplied by the wavelength since it's 1/Angstrom?
            # flux_fact*=wvlength 
            # or the bandwidth ?
            flux_fact*= bandwidth*u.erg/(u.second*u.cm*u.cm)
            # to convert into erg cm–2 sec–1 Hz–1
            flux_fact = flux_fact.to("(erg*Hz)/ cm2")#.value
            ZP_mag    = get_header(path_llm,"ZP") #AB mag 
        else:
            #pragma no cover -> in these cases, the radius we are considering is at the image pos 
            # meaning that the light enclosed would be "spoiled" by the qso light
            raise RuntimeError("still to implement")
            """
            # from HST images 
            lens_light_pix = lens_light(sett,unconvolved=True,kwres=kw_res,bandmodel=None)
            path_data = setting.data_path+"/"+setting.image_name
            if "475" in get_filter(sett) or "814" in get_filter(sett):
                # UVIS: https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/uvis-encircled-energy
                # Flux = Fr*photoflam/EE(r)
                PHOTOFLAM = get_header(path_data,"PHTFLAM1") # Jy * sec/count 
                EE_r = get_EE(sett,radius)
            #ZP_mag    = 8.9 # from https://www.stsci.edu/documents/dhb/web/c20_nicdataanal.fm2.html
            photonu   = get_header(path_data,"PHOTFNU")*u.second/1 # Jy * sec/count
            """
        #flux  = flux_in_grid(lens_light_pix,grid)*flux_fact
        flux  = flux_in_grid(lens_light_pix,grid)
        mag   = -2.5*np.log10(flux.value) + ZP_mag 
        print("TESTING: Mag whithin radius ", format(radius,".2f"),"\": ",format(mag,".2f")," mag (AB) calc w ZP")
        flux  = flux_in_grid(lens_light_pix,grid)*flux_fact
        mag   = -2.5*np.log10(flux.value) 
        print("TESTING: Mag whithin radius ", format(radius,".2f"),"\": ",format(mag,".2f")," mag (AB) calc w conversion in energy/sec")
        
        print("Mag whithin radius ", format(radius,".2f"),"\": ",format(mag,".2f")," mag (AB)")
        
        luminosity = 4*np.pi*(cosmoLd**2)*flux # luminosity distance
        luminosity = luminosity.to("solLum")
        print("Luminosity within radius ", format(radius,".2f"),"\": ",format(luminosity.value,".2E")," L_sun")
        
        Mass_to_light_ratio = mass_discreet/luminosity
        print("Mass-to-light ratio within radius ",format(radius,".2f"),"\": ",format(Mass_to_light_ratio.value,".2E")," M_sun/L_sun")
    success(sys.argv[0])
