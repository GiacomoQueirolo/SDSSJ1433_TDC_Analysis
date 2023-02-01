from tools import *
from get_res import *
from image_manipulation import *
from pos_img import create_pos_img
from create_fits import lens_light
from input_data import init_lens_model
from conversion import conv_radec_to_xy,conv_xy_to_radec

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import default_cosmology
from lenstronomy.LensModel.lens_model import LensModel
    
class Grid_Class():
    # create squared grid with of centers of pixels
    # given center, "radius" (bc it can be made circular) 
    def __init__(self,center_x,center_y,radius,precision=0.1,type="pixel",setting=None):
        self.type       = type # "pixel" or "radec" 
        self.center     = [center_x,center_y]
        self.radius     = radius
        self.precision  = precision
        self.grid       = self.create_grid()
        self.grid_area  = self.get_grid_area()
        self.pixel_grid = None
        self.radec_grid = None
        if setting:
            self.setting = get_setting_module(setting).setting()
        if self.type=="pixel":
            self.pixel_grid = self.grid
        else:
            self.pixel_grid = self.get_transf_grid("to_pixel")
        
        if self.type=="radec":
            self.radec_grid = self.get_transf_grid("to_arcsec") 
        else:
            self.radec_grid = self.grid
    
    def create_grid(self):
        n  = int(2*self.radius/self.precision)
        n  = max([1,n]) 
        center_x,center_y = self.center
        dx = np.linspace(center_x-self.radius,center_x+self.radius,n)
        dy = np.linspace(center_y-self.radius,center_y+self.radius,n)
        grid_centers = [dx,dy]
        grid_of_points = np.transpose(np.meshgrid(*grid_centers))
        return grid_of_points

    def get_grid_area(self):
        # easy way out: assume constant area:
        # grid ~ = [[[a,b],[a,c]],[[d,b],[d,c]]] for 4 points
        y_values = self.grid[0].T[1]  
        x_values = self.grid.T[0][0] 
        dy = [ y_values[i]-y_values[i-1] for i in range(1,len(y_values))]
        dx = [ x_values[i]-x_values[i-1] for i in range(1,len(x_values))]
        prec = 1e-5
        if np.all(np.abs(dx-dx[0])<prec) and np.all(np.abs(dy-dy[0])<prec):
            area = dx[0]*dy[0]
            self.grid_area = np.ones(np.shape(self.grid)[:-1])*area
            return self.grid_area
        else:
            print(dx,dy)
            # pragma: no cover
            raise 
    
    def _d_center(self,grid_flat):
        #d_i = np.sqrt((grid_flat[i][0]-self.center[0])**2+(grid_flat[i][1]-self.center[1])**2)
        #return d_i
        for gi in grid_flat:
            yield np.sqrt((gi[0]-self.center[0])**2+(gi[1]-self.center[1])**2)
        
    def circularise_flat_grid_with_area(self):
        grid_flat = self.grid.reshape(len(self.grid)**2,2)
        grid_area_flat = self.grid_area.reshape(len(self.grid_area)**2)
        """grid = []
        area = []
        for i in range(len(grid_flat)):
            d_i_center = np.sqrt((grid_flat[i][0]-self.center[0])**2+(grid_flat[i][1]-self.center[1])**2)
            if d_i_center<=self.radius:
                grid.append(np.array(grid_flat[i]).tolist())
                area.append(np.array(grid_area_flat[i]).tolist())
        ###
        grid  = np.array([grid_flat[i] for i in range(len(grid_flat)) if self._d_center(grid_flat,i)<=self.radius])
        area  = np.array([grid_area_flat[i] for i in range(len(grid_area_flat)) if self._d_center(grid_flat,i)<=self.radius])
        """
        d_center = self._d_center(grid_flat) 
        select_pix = [True  if d<=self.radius else False for d in d_center ]
        grid = np.array(grid_flat[select_pix])
        area = np.array(grid_area_flat[select_pix])
        #area = np.array([ for i,d in zip(range(len(grid_area_flat),d_center) if self._d_center(grid_flat,i)<=self.radius])
        self.circ_flat      = grid
        self.area_circ_flat = area
        return self.circ_flat,self.area_circ_flat
    

    def get_transf_grid(self,transf="to_pixel"):
        if transf=="to_pixel" and hasattr(self,"pixel_grid"):
            return self.pixel_grid
        elif transf=="to_arcsec" and hasattr(self,"radec_grid"):
            return self.radec_grid
        grid_trnsf = np.ones_like(self.grid)
        """for i in range(len(self.grid)):
            for j in range(len(self.grid[i])):
                if   transf=="to_pixel":
                    grid_trnsf[i][j] = np.reshape(conv_radec_to_xy(self.setting,*self.grid[i][j]),2)
                elif transf=="to_radec":
                    grid_trnsf[i][j] = np.reshape(conv_xy_to_radec(self.setting,*self.grid[i][j]),2)
        """
        if transf not in ["to_pixel","to_arcsec"]:
            raise RuntimeError("trasnformation method cannot be "+str(transf)+" but either to_pixel or to_arcsec")
        grid_trnsf = [[np.reshape(conv_radec_to_xy(self.setting,*gij),2) if transf=="to_pixel" else np.reshape(conv_xy_to_radec(self.setting,*gij),2) for gij in gi] for gi in self.grid]
        return grid_trnsf
    
    def convert_grid(self):
        pixtoarcsec = get_pixscale(self.setting)
        if self.type=="radec":
            #convert to pixel
            self.type      = "pixel"
            self.center    = np.reshape(conv_radec_to_xy(self.setting,*self.center),2).tolist()#conv_radec_to_xy(sets,ra=max_ra,dec=max_dec)f
            self.radius    = self.radius/pixtoarcsec
            self.precision = self.precision/pixtoarcsec
            self.grid      = self.get_transf_grid("to_pixel")
            self.grid_area = self.get_grid_area()
            if hasattr(self,"circ_flat"): 
                self.circularise_flat_grid_with_area()
            return 0
        elif self.type=="pixel":
            #convert to radec
            self.type      = "radec"
            self.center    = np.reshape(conv_xy_to_radec(self.setting,*self.center),2).tolist()
            self.radius    = self.radius*pixtoarcsec
            self.precision = self.precision*pixtoarcsec
            self.grid      = self.get_transf_grid("to_arcsec")
            self.grid_area = self.get_grid_area()
            if hasattr(self,"circ_flat"): 
                self.circularise_flat_grid_with_area()
            return 0
        
    def degrade_to_pixel(self):
        if self.type=="radec":
            self.convert_grid()
        grid_int = np.rint(self.grid).astype(int)
        # we eliminate the multiple pixels
        # we also loose the order of the pixel, but it's not important
        flat_grid_int = np.array(grid_int).reshape(len(grid_int)**2,2)
        set_grid = set(tuple(gp) for gp in flat_grid_int)
        # it can't be not flat as it has deleted some points
        self.flat_grid_degraded = np.array(list(set_grid))
        # the area of a pixel is 1 pixel in pixel units :
        self.flat_area_degraded = np.ones_like(flat_grid_degraded) 
        return self.flat_grid_degraded
            
    def plot_grid(self,color="green",ax=None,circ=False):
        if not ax:
            fig,ax=plt.subplots()
        pix_width  = self.precision
        pix_height = pix_width
        if not circ:
            flat_grid   = self.grid.reshape(len(self.grid)**2,2) 
        elif hasattr(self,"circ_flat"):
            flat_grid   = self.circ_flat
        else:
            flat_grid,_ = self.circularise_flat_grid_with_area() 
        for i in range(len(flat_grid)):
            anchor_p = [flat_grid[i][0]-pix_width/2.,flat_grid[i][1]-pix_height/2.]
            rect = patches.Rectangle(anchor_p, pix_width - 0.01*self.precision, pix_height - 0.01*self.precision,color = color)
            ax.add_patch(rect)
        return ax
    
        


def get_cosmo_prm(setting):
    z_source  = setting.z_source
    z_lens    = setting.z_lens

    cosmo     = default_cosmology.get()
    cosmo_dd  = cosmo.angular_diameter_distance(z_lens)
    cosmo_ds  = cosmo.angular_diameter_distance(z_source)
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source)
    
    Sigma_crit = const.c*const.c*cosmo_ds/(4*np.pi*const.G*cosmo_dds*cosmo_dd)
    Sigma_crit = Sigma_crit.to("Msun/(kpc kpc)")
    return Sigma_crit,cosmo_dd



def get_kappa_discreet(setting,radius,lens_model,kwlens,dr=0.1):
    center_lens_x,center_lens_y = kwlens[0]["center_x"],kwlens[0]["center_y"] 
    grid      = Grid_Class(center_x=center_lens_x,center_y=center_lens_y,radius=radius,precision=dr,type="radec",setting=setting)
    grid_circ_flat,grid_area_circ_flat = grid.circularise_flat_grid_with_area()
    if not (np.shape( grid_circ_flat)[:-1]==np.shape(grid_area_circ_flat)):
        raise
    kappa_flat  = lensModel.kappa(*grid_circ_flat.T,kwargs=kwlens)
    return kappa_flat,grid
    


def get_mass_discreet(setting,radius,lens_model,kwlens,Sigma_crit,cosmo_dd,dr=0.1,ret_grid=False):
    # dr= precision in arcsec
    # Mass(r) = 2 pi Sigma_crit * integr_0^r dr kappa(r) * r
    # the integrand is then:
    # or rather, since it's discretised anyway:
    # M(r) = Sigma_crit * Sum_k kappa(k)*Area(k) with k indicating the grid cells within the radius
    kappa_discreet, grid =  get_kappa_discreet(setting,radius,lens_model,kwlens,dr=dr) 
    area_grid      = grid.area_circ_flat
    area_grid_Mpc2 = area_grid*cosmo_dd*cosmo_dd*u.arcsec.to("rad")*u.arcsec.to("rad")
    Mass = np.sum(kappa_discreet*area_grid_Mpc2)*Sigma_crit 
    Mass = Mass.to("Msun")
    #mass_3D   = lensModel.lens_model.mass_3d(r=radius,kwargs=kwlens["kwargs_lens"],bool_list=[0]) #bc shear do not have the mass and SIS is ignored (see note may 23th '22)
    #mass_fact = u.arcsec.to("rad") ** 2 * cosmo_dd * cosmo_ds / cosmo_dds * const.c ** 2 / (4 * np.pi * const.G)

    #mass = mass_3D*mass_fact
    #mass = mass.to("Msun") #solar masses
    return Mass,grid 

#################################
# get light from pixellated map #
#################################

def flux_in_grid(lens_light_pix,grid):
    if grid.type=="radec":
        grid = grid.convert_grid()
    # now the problem is, the grid is not discretised in pixels
    flat_grid_degraded = grid.degrade_to_pixel()
    # we now have a grid with the coord of the pixels
    flux = np.sum([lens_light_pix[i[0]][i[1]] for i in flat_grid_degraded])
    return flux
    
if __name__=="__main__":
    present_program(sys.argv[0])

    parser = argparse.ArgumentParser(description="Mass, Light and Mass-to-Light ratio within a certain radius for given models")
    #parser.add_argument("-PR", "--plot_wrt_radius",dest="plot_wrt_radius", default=False,action="store_true",
    #                    help="Plot with respect to radius")
    parser.add_argument("-r", "--radius", type=int, dest="radius", default=-1,
                        help="Radius at which to calculate the mass. Default: approx max radius of the lensed image position")                 
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting files to consider")
    
    args      = parser.parse_args()
    #PR        = args.plot_wrt_radius
    radius    = args.radius
    setting_names  = args.SETTING_FILES
    ###
    
    dr = 0.001
    
    ###
    to_redo = False

    for sett in setting_names:
        print_setting(sett)
        setting = get_setting_module(sett).setting()
        kw_res  = get_kwres(setting)["kwargs_results"]
        kwlens  = kw_res["kwargs_lens"]
        lensModel = init_lens_model(setting)
        SigmaCrit,cosmoDd = get_cosmo_prm(setting)
        
        print("dr :",dr*cosmoDd.to("kpc")*u.arcsec.to("rad") )
        if radius==-1 or to_redo:
            print("Consider the mass enclosed in a radius of approximately at the image positions, centered on the center of the main lens")
            #print("from lenstronomy: \"This function assumes spherical symmetry/ignoring the eccentricity.\" ")
            ps_res = kw_res["kwargs_ps"][0]
            ra_image, dec_image = ps_res["ra_image"],ps_res["dec_image"] # A is ~0,0
            darcsec = []
            for i in range(len(ra_image)):
                darcsec.append(max([ np.sqrt( (ra_image[i]-x)**2 + (dec_image[i]-y)**2) for x,y in zip(ra_image,dec_image)] ))
            radius  = max(darcsec)/2.
            to_redo = True
        mass_discreet,grid = get_mass_discreet(setting,radius,lensModel,kwlens,SigmaCrit,cosmoDd,dr=dr,ret_grid=True)
        print("Mass enclosed in radius ", format(radius,".2f"),"\": ",format(mass_discreet.value,".2E")," Msun")
        # plot pixellated grid: 
        # TO TEST
        ax_im = create_pos_img(setting,RP=False,RM=False,save=False,also_fits=False)
        if grid.type=="radec":
            grid.convert_grid()
        ax_im_pix = grid.plot_grid(ax=ax_im,circ=True)
        plt.savefig(get_savefigpath(setting)+"/grid_plot.png")
        
        flux_unit    = 1/u.second #count per sec
        #https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints
        #https://hst-docs.stsci.edu/wfc3dhb/chapter-9-wfc3-data-analysis/9-1-photometry
        if check_SUB(sett): # this is SUB=> model from isophotes
            lens_light_pix = get_lens_light_model(sett)
            path_llm = setting.data_path+"/"+setting.lens_light_model_name
            #ZP_mag   = get_header(path_llm,"ZP")
            #photonu  = get_header(path_llm,"PHOTFNU")*u.second/1 # ~Jy~ * sec/count -> not Jy ? must somehow already be taken into account  
            #if "475" in get_filter(sett) or "814" in get_filter(sett):
            #   UVIS : https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/uvis-encircled-energy
            #   IR : https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/ir-photometric-calibration
            #   for both Flux = Fr*photoflam/EE(r)
            wvlenght  = get_header(path_llm,"PHOTPLAM")# Angstrom
            bandwidth = get_header(path_llm,"PHOTBW")# Angstrom
            PHOTOFLAM = get_PHOTOFLAM(setting)
            EE_r      = get_EE(setting,radius)
            flux_fact = PHOTOFLAM/EE_r # ergs/cm2/Angs/e-
            # should this be multiplied by the wavelenght since it's 1/Angstrom?
            # flux_fact*=wvlenght 
            # or the bandwidth ?
            flux_fact*=bandwidth*u.erg/(u.second*u.cm2)
            ZP_mag    = get_header(path_llm,"ZP")
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
        flux  = flux_in_grid(lens_light_pix,grid)*flux_fact
        mag   = -2.5*log10(flux) + ZP_mag 
        print("Mag whithin radius ", format(radius,".2f"),"\": ",format(mag,".2f")," mag (AB)")
        
        cosmo = default_cosmology.get()
        dist  = cosmo.luminosity_distance(setting.z_lens) 
        luminosity = 4*np.pi*(dist**2)*flux
        luminosity = luminosity.to("solLum")
        print("Luminosity within radius ", format(radius,".2f"),"\": ",format(luminosity.value,".2f")," L_sun")
        
        Mass_to_light_ratio = mass_discreet/luminosity
        print("Mass-to-light ratio within radius ",format(radius,".2f"),"\": ",format(Mass_to_light_ratio.value,".2f")," M_sun/L_sun")
    success(sys.argv[0])
