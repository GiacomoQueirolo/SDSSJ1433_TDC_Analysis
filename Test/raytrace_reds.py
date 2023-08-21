# following suggestion of Sluse,
# I try raytracing back the position of the "red sources"
# and back again into the image plane
# using lens model given by either the f140 or f160 
# and see where it land/if it's compatible with the other weak
# signal at the bottom
import sys
import argparse
import numpy as np
from copy import copy
from Utils.tools import *
from Utils.get_res import * 
from Data.input_data import *
from Data.reg_mask import str_reg
from MassToLight.grid_class import dist
from Data.conversion import conv_xy_to_radec,conv_radec_to_xy

from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from astropy.cosmology import FlatLambdaCDM

print("To update, look up raytrace_reds.py in usm machines")
#exit()
# from f160w 
# approx pix pos X,Y from DS9 (corrected automatically)
pos_R1_px = np.array([62.,64.98])
pos_R2_px = np.array([39.875,42.17])
pos_R1 = conv_xy_to_radec(x=pos_R1_px[0],y=pos_R1_px[1],setting="f160w_CPRP_PLLI")
pos_R2 = conv_xy_to_radec(x=pos_R2_px[0],y=pos_R2_px[1],setting="f160w_CPRP_PLLI")
# rename the object as Red1 and Red2
def get_xy_images(lensModel,kw_lens,pos_R=pos_R1):
    pert_main_source_plane = lensModel.ray_shooting(*pos_R[0],*pos_R[1],kw_lens)
    LensEqSolver = LensEquationSolver(lensModel)
    new_x,new_y  = LensEqSolver.image_position_from_source(*pert_main_source_plane,kwargs_lens=kw_lens)
    return new_x,new_y
    

def get_xy_reg(settings,lensModel,kw_lens,pos_R=pos_R1,mag_rad=True):
    new_x,new_y = get_xy_images(lensModel,kw_lens,pos_R=pos_R)
    rad = 4
    mag = 1
    str_reg_pert = copy(str_reg) 
    if type(settings) is list:
        str_reg_pert_list = []
        for setting in settings:
            for ra,dec in zip(new_x,new_y):
                x,y = conv_radec_to_xy(setting=setting,ra=ra,dec=dec)
                if mag_rad:
                    mag = lensModel.magnification(ra,dec,kw_lens)
                str_reg_pert+=f"circle({x[0]+1},{y[0]+1},{rad*mag}) # color=green dash=1\n"
            str_reg_pert_list.append(str_reg_pert)
        return str_reg_pert_list
    else:
        for ra,dec in zip(new_x,new_y):
            if mag_rad:
                mag = lensModel.magnification(x,y,kw_lens)
            x,y = conv_radec_to_xy(setting=setting,ra=ra,dec=dec)
            str_reg_pert+=f"circle({x[0]+1},{y[0]+1},{rad*mag}) # color=green dash=1\n"
        return str_reg_pert

            
# further test: see if the positions can be better fitted by changing redshift

def get_dist_R(new_pos_Ri,pos_Ri=[pos_R1,pos_R2]):
    dist_R =  [min([dist(new_pos_Ri_i,pos_Ri_i) for new_pos_Ri_i in new_pos_Ri]) for pos_Ri_i in pos_Ri]
    return dist_R

def mod_w_z(kw_lens,z_lens,z_source,z_source_prime,cosmo=None):
    # theta_E and shear strenght are changed by a factor when the redshift changes
    # theta_E \propto sqrt(D_LS/D_L*D_S) -> theta_E(z') = theta_E*sqrt(D_LS(z')/*D_S(z'))*sqrt(D_S(z)/*D_LS(z))
    #fact :sqrt(D_LS(z')/*D_S(z'))*sqrt(D_S(z)/*D_LS(z))
    if not cosmo:
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # cosmo from https://academic.oup.com/mnras/article/474/3/3391/4644836, Agnello 2017
    cosmo_ds  = cosmo.angular_diameter_distance(z_source)
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source)
    new_cosmo_ds  = cosmo.angular_diameter_distance(z_source_prime)
    new_cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source_prime)
    fact = np.sqrt(cosmo_dds*new_cosmo_ds/(cosmo_ds*new_cosmo_dds))
    kw_lens[0]["theta_E"]   = kw_lens[0]["theta_E"]*fact
    kw_lens[1]["theta_E"]   = kw_lens[1]["theta_E"]*fact
    kw_lens[2]["gamma_ext"] = kw_lens[2]["gamma_ext"]*fact  
    return kw_lens
    

def fit_redshift(setting,kw_lens,range_zpert=None,lens_model_list=None,pos_R=pos_R1,nzp=10):
    range_zpert[0] = min([setting.z_lens,*range_zpert])
    if kw_lens is None:
        kw_lens = get_kwres(setting_name=lens_model_setting)["kwargs_results"]["kwargs_lens"]
    if lens_model_list is None:
        lens_model_list = init_lens_model_list(setting)
    dist_R = []
    #zls    = np.linspace(*range_zlens,nzl)
    zps    = np.linspace(*range_zpert,nzp)
    
    lensModel = LensModel(lens_model_list=lens_model_list)
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3) 
    for zp in zps:
        print("z_pert",zp)#, z_lens=zl, z_source=setting.z_source)
        new_kw_lens   = mod_w_z(kw_lens,setting.z_lens,setting.z_source,z_source_prime=zp,cosmo=cosmo)
        new_x,new_y   = get_xy_images(lensModel=lensModel,kw_lens=new_kw_lens,pos_R=pos_R)

        dist1,dist2 = get_dist_R(np.transpose([new_x,new_y]))
        dist_R.append([dist1,dist2]) #dist1,dist2
    dist_R = np.transpose(dist_R) # shape: z,2
    for d in dist_R:
        plt.scatter(zps,d[0],label="d(R1)")
        print(d)
        if min(dist_R.T[0])==d[0]:
            plt.axvline(d,label="min d(R1)")
        if min(dist_R.T[1])==d[1]:
            plt.axvline(d,label="min d(R2)")
        plt.scatter(zps,d[1],label="d(R2)")
        plt.scatter(zps,np.sum(d),label="d(R1)+d(R2)")
    plt.legend()
    filename = get_savefigpath(setting)+"/err_vs_redshift_R12.pdf"
    plt.savefig(filename)
    print(f"Saved file {filename}")

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(prog=sys.argv[0],description="Ray-trace perturber",formatter_class=CustomFormatter)
    parser.add_argument('-lms','--lens_model_setting',type=str,dest="lens_model_setting",default="f140w_CPRP_PLLI",help="Setting from which consider the lens model")
    parser.add_argument('-rms','--red_model_setting',type=str,dest="red_model_setting",default="f160w_CPRP_PLLI",help="Setting where the 'red sources' are observed (f160w)")
    args = parser.parse_args()
    lens_model_setting = get_setting_module(args.lens_model_setting,1)
    red_model_setting  = get_setting_module(args.red_model_setting,1)
    settings           = [lens_model_setting,red_model_setting]
    lensModel    = init_lens_model(lens_model_setting)
    kw_lens      = get_kwres(setting_name=lens_model_setting)["kwargs_results"]["kwargs_lens"]
    str_reg_list = get_xy_reg(settings=settings,lensModel=lensModel,kw_lens=kw_lens,pos_R=pos_R1)
    for str_reg,sett in zip(str_reg_list,settings):
        reg_file     = get_savefigpath(setting=sett)+"/raytraced_pos_red.reg"
        with open(reg_file, "w") as text_file:
            print(str_reg, file=text_file)
        print(f"Saved file {reg_file}")
    fit_redshift(lens_model_setting,kw_lens)
    success(sys.argv[0])

