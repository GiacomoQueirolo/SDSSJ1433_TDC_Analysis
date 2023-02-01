from tools import *
from conversion import conv_radec_to_xy,conv_xy_to_radec
from image_manipulation import get_rotangle,fits_with_copied_hdr, get_numPix

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
#import timeit
 

def dist(xy0,xy1):
    return np.sqrt( (xy0[0]-xy1[0])**2 + (xy0[1]-xy1[1])**2 ) 
    

class Point():
    def __init__(self,x,y,center,precision,setting,input_type="radec"):
        self.setting = get_setting_module(setting,sett=True)
        try:
            self.center_radec = center.ra,center.dec
            self.center_px    = center.x_pix,center.y_pix
        except:        
            if input_type=="radec":
                self.center_radec = center
                self.center_px    = np.reshape(conv_radec_to_xy(self.setting,*center),2).tolist()
            elif input_type=="pixel":
                self.center_px    = center
                self.center_radec = np.reshape(conv_xy_to_radec(self.setting,*center),2).tolist()
                
        if input_type=="radec":
            self.precision_radec  = precision
            self.precision_px     = precision/self.setting.pix_scale
            self.ra,self.dec      = x,y
            self.x_pix,self.y_pix = np.reshape(conv_radec_to_xy(self.setting,self.ra,self.dec),2).tolist()
        elif input_type=="pixel":
            self.precision_px     = precision
            self.precision_radec  = precision*self.setting.pix_scale
            self.x_pix,self.y_pix = x,y
            self.ra,self.dec      = np.reshape(conv_xy_to_radec(self.setting,self.x_pix,self.y_pix),2).tolist()

        self.ra,self.dec      = self.get_radec()
        self.x_pix,self.y_pix = self.get_pixel()
        self.dist_from_center_pix   = self.get_dist_from_center("pixel")
        self.dist_from_center_radec = self.get_dist_from_center("radec")
        
    def get_pixel(self):
        return self.x_pix,self.y_pix
        
    def get_pixel_int(self):
        if hasattr(self,"x_pix_int") and hasattr(self,"y_pix_int"):
            return self.x_pix_int,self.y_pix_int
        x_pix,y_pix = self.get_pixel()
        self.x_pix_int,self.y_pix_int = np.round(x_pix,0).astype(int),np.round(y_pix,0).astype(int) 
        return self.x_pix_int,self.y_pix_int
        
    def get_radec(self):
        return self.ra,self.dec
        
    def get_dist_from_center(self,output_type="pixel"):
        if output_type not in ["pixel","radec"]:
            raise ValueError("Output either in pixel or radec, not "+str(output_type))
        if output_type=="pixel":
            if not hasattr(self,"dist_from_center_pix"):
                self.dist_from_center_pix = dist(self.center_px,self.get_pixel())
            return  self.dist_from_center_pix
        elif output_type=="radec":
            if not hasattr(self,"dist_from_center_radec"):
                self.dist_from_center_radec = dist(self.center_radec,self.get_radec())
            return  self.dist_from_center_radec
    def dist_from(self,P2):
        return dist_btw_points(self,P2)
    
    def __add__(self,shift):
        kwarg_point = {"center":self.center_radec, "precision":self.precision_radec, 
                     "setting":self.setting, "input_type":"radec"}
        if not isinstance(shift,Point) and type(shift,list):
            shift = Point(*shift,**kwarg_point)
            print("Warning: Assuming shift in radec")
        elif shift.center_radec!=self.center_radec:
            shift_center = [self.center_ra+shift.center_ra,self.center_dec+shift.center_dec]
            shift_ra  = shift.ra  + shift_center[0]
            shift_dec = shift.dec + shift_center[1]
            shift     = Point(shift_ra,shift_dec,**kwarg_point)
        ra  = self.ra+shift.ra
        dec = self.dec+shift.dec
        return Point(ra,dec,**kwarg_point)

        
"""
    # to consider overwriting "P1+P2"
    # but complicated if they don't have the same center
    def __add__(self, P2):
        if self.center_radec ==P2.center_radec:
            shift_radec= [0,0]
        el
        x = self.ra  + P2.ra
        y = self.dec + P2.dec
        return Point(x)
""" 
class Square(Point):
    def __init__(self,x,y,center,precision,setting,input_type="radec",radius=None):
        self.Point      = Point(x,y,center,precision,setting,input_type)
        super().__init__(x,y,center,precision,setting,input_type)
        self.center     = center
        self.area_px    = self.precision_px*self.precision_px
        self.area_radec = self.precision_radec*self.precision_radec
        if radius:
            if not isinstance(radius,Length):
                radius = Length(radius,precision,setting,input_type)
            cent_dist = self.dist_from_center()
            if cent_dist.length_pix<radius.length_pix:
                self.inside_rad = True
            else:
                self.inside_rad = False
                
    def dist_from_center(self):
        dist_pix = self.get_dist_from_center("pixel") # dist([self.x_pix,self.y_pix],[self.center.x_pix,self.center.y_pix])
        return Length(dist_pix,self.setting,input_type="pixel")
    
    def __add__(self,shift):
        print("Warning: we are assuming this is a shift and moreover a shift of the whole grid")
        new_Point = self.Point+shift
        new_Square = Square(new_Point.ra,new_Point.dec,new_Point.center_radec,new_Point.precision_radec,
                            new_Point.setting,input_type="radec")
        # we assume is a shift of the whole grid: it doesn't change wheter this is within or outside the radius
        if hasattr(self,"inside_rad"):
            new_Square.inside_rad = self.inside_rad
        return new_Square
            
            

class Length():
    def __init__(self,length,setting,input_type="radec"): 
        self.setting = get_setting_module(setting,sett=True)
        if input_type == "radec":
            self.length_radec = length 
            self.length_pix   = length/self.setting.pix_scale
        elif input_type=="pixel":
            self.length_radec = length*self.setting.pix_scale
            self.length_pix   = length
        #self.length_pix   = length/precision
        #_ra,_dec  =  np.reshape(conv_xy_to_radec(self.setting,0,self.length_pix),2).tolist()
        #self.length_radec = dist([_ra,_dec],[0,0])
class Shift(Length):
    def __init__(self,dx,dy,setting,input_type="radec"):
        if input_type=="radec":
            self.dra,self.ddec      = dx,dy
            self.dx_pix,self.dy_pix = np.reshape(conv_radec_to_xy(self.setting,self.dra,self.ddec),2).tolist()
        elif input_type=="pixel":
            self.dx_pix,self.dy_pix = dx,dy
            self.dra,self.ddec      = np.reshape(conv_xy_to_radec(self.setting,self.dx_pix,self.dy_pix),2).tolist()
        length=np.sqrt(dx**2+dy**2)
        super().__init__(length,setting,input_type)
        

def dist_btw_points(P1,P2):
    if not isinstance(P1,Point) or not isinstance(P2,Point):
        raise ValueError("The two imputs have to be instances of class Point. Else use dist(xy0,xy1)")
    if get_setting_name(P1.setting)!=get_setting_name(P2.setting):
        raise ValueError("The points given must have the same setting")
    ra1,dec1 = P1.get_radec()
    ra2,dec2 = P2.get_radec()
    dist_radec = dist([ra1,dec1],[ra2,dec2])
    return Length(dist_radec,P1.setting,input_type="radec")
    
class Grid_Class():
    # create squared grid with centers of pixels
    # given center, "radius" (bc it can be made circular) 
    def __init__(self,center,radius,precision,setting,input_type="radec",r_mask_cnt=None):
        self.Center = Point(*center,[0,0],precision,setting,input_type)
        self.Radius = Length(radius,setting,input_type)
        if r_mask_cnt:
            self.Rmask_cnt   = Length(r_mask_cnt,setting,input_type)
        else:
            self.Rmask_cnt   = None
        self.precision_px    = self.Center.precision_px
        self.precision_radec = self.Center.precision_radec
        self.setting         = get_setting_module(setting,sett=True)
        self._sqrt_grid_radec = self._create_squared_grid_of_radec_points()
        self._sqrt_grid       = self._get_grid()
        self.circ_grid_flat   = self.get_circ_grid_flat()
        
    def _create_squared_grid_of_radec_points(self):
        # let's first consider it in radec
        rd_radec = self.Radius.length_radec
        n  = int(2*rd_radec/self.precision_radec)
        n  = max([1,n]) 
        center_ra,center_dec = self.Center.ra,self.Center.dec
        
        dra  = np.linspace(center_ra-rd_radec,center_ra+rd_radec,n)
        ddec = np.linspace(center_dec-rd_radec,center_dec+rd_radec,n)
        grid_centers_radec   = [dra,ddec]
        grid_of_points_radec = np.transpose(np.meshgrid(*grid_centers_radec))
        return grid_of_points_radec
        
    def _get_grid(self):
        if hasattr(self,"_sqrt_grid"):
            return self._sqrt_grid
        kwargs = {"center":self.Center,"precision":self.precision_radec,"setting":self.setting,\
                  "input_type":"radec","radius":self.Radius}
        grid = [[Square(*gij,**kwargs) for gij in gi] for gi in self._sqrt_grid_radec]
        return grid
    
    def _flatten(self,grid):
        try:
            return np.array(grid).reshape(len(grid)**2,2)
        except ValueError: # 
            return np.array(grid).flatten()   
        
    def get_grid_flat(self):
        return self._flatten(self._get_grid())
    
    def get_circ_grid_flat(self):
        if hasattr(self,"circ_grid_flat"):
            return self.circ_grid_flat
        grid_select = [[gij.inside_rad for gij in gi] for gi in self._sqrt_grid]
        grid_circ = self._flatten(self._sqrt_grid)[self._flatten(grid_select).tolist()]
        if hasattr(self,"Rmask_cnt"):
            grid_circ = self._mask_center(grid_circ,self.Rmask_cnt)
        return np.array(grid_circ)
        
    def _get_degraded_pixel_grid_flat(self):
        """
        def _get_degraded_pixel_grid_flat(self,flat_grid=None):
        flat_grid_int = [cgf.get_pixel_int() for cgf in flat_grid] 
        # we eliminate the multiple pixels
        # we also loose the order of the pixel, but it's not important
        #flat_grid_int = grid_pix_int
        set_grid = set(tuple(gp) for gp in flat_grid_int)
        # it can't be not flat as it has deleted some points
        # the area of a pixel is 1 pixel in pixel units :
        #flat_area_degraded = np.ones_like(flat_grid_degraded) 
        return np.array(list(set_grid)) # those are only points not instances of Points
        
        """
        # other idea: return all pixels within the radius
        # !: not 1to1 relation btw the pixels and radecs
        Xc,Yc = self.Center.get_pixel()
        Rc    = self.Radius.length_pix
        topright_grid   = np.round([Xc+Rc,Yc+Rc],0).astype(int)
        bottomleft_grid = np.round([Xc-Rc,Yc-Rc],0).astype(int)
        
        xx = np.arange(bottomleft_grid[0],topright_grid[0])
        yy = np.arange(bottomleft_grid[1],topright_grid[1])
        
        grid_full = np.transpose(np.meshgrid(xx,yy))
        grid_full = grid_full.reshape(grid_full.shape[0]**2,2)
        grid_cut  = []
        for xyi in grid_full:
            if dist(xyi,[Xc,Yc])<=Rc:
                if masked=True:
                    if dist(xyi,)
                grid_cut.append(xyi.tolist())
        return np.array(grid_cut)
    
    def get_degraded_pixel_grid_flat(self):
        if hasattr(self,"degraded_pixel_grid_flat"):
            return self.degraded_pixel_grid_flat
        # self.degraded_pixel_grid_flat = self._get_degraded_pixel_grid_flat(self.circ_grid_flat)
        self.degraded_pixel_grid_flat = self._get_degraded_pixel_grid_flat()
        return self.degraded_pixel_grid_flat 
    
    def get_area_square(self,output_type="radec"):
        if output_type not in ["pixel","radec"]:
            raise ValueError("Output either in pixel or radec, not "+str(output_type))
        # they are the same for every square by construction
        if output_type=="pixel":
            return self._sqrt_grid[0][0].area_px
        elif output_type=="radec":
            return self._sqrt_grid[0][0].area_radec
            
    def get_area_grid(self,output_type="radec"):
        if output_type not in ["pixel","radec"]:
            raise ValueError("Output either in pixel or radec, not "+str(output_type))
        # still the same area for each square
        return self.get_area_square(output_type)*len(self.get_circ_grid_flat())
    
    def plot_grid(self,color="green",ax=None,grid_to_plot="circ_grid_flat",alpha=.7):
        if not ax:
            fig,ax=plt.subplots()
        if grid_to_plot not in ["sqrt_grid","circ_grid_flat","degraded_pixel_grid_flat"]:
            raise RuntimeError
        else:
            if grid_to_plot =="sqrt_grid":
                class_point=True
                flat_grid = self.get_grid_flat()
                pix_width = self.precision_px
                angle     = get_rotangle(self.setting)
            elif grid_to_plot=="circ_grid_flat":
                class_point=True
                flat_grid = self.get_circ_grid_flat()
                pix_width = self.precision_px
                angle     = get_rotangle(self.setting)
            
            elif grid_to_plot=="degraded_pixel_grid_flat":
                class_point=False # those are just coordinates
                pix_width = 1
                angle     = 0
        pix_height = pix_width
        PC = []
        for i in range(len(flat_grid)):
            if not class_point:
                anchor_p = [flat_grid[i][0]-pix_width/2.,flat_grid[i][1]-pix_height/2.]
            else:
                anchor_p = [flat_grid[i].x_pix-pix_width/2.,flat_grid[i].y_pix-pix_height/2.]
            rect = patches.Rectangle(anchor_p, pix_width - 0.05*pix_width, pix_height - 0.05*pix_height,color = color,angle=angle)
            PC.append(rect)
        ax.add_collection(PatchCollection(PC,facecolor=color,alpha=alpha,edgecolor=None))
        # remember to do "ax.plot()" in order to plot it!
        return ax 
    
    def produce_degraded_fits(self,name="degraded_grid",verbose=True):
        # in order to be sure that the problem with the grid is not due to plotting
        # i now produce the fits file to compare w the original image
        # for the degraded pixel grid
        grid = self.get_degraded_pixel_grid_flat()
        main_image_path = self.setting.data_path+"/"+self.setting.image_name
        numpix=get_numPix(self.setting,twodim=True)
        mask = np.zeros(numpix)
        for px in grid:
            mask[px[0]][px[1]] = 1
        svpath = self.setting.data_path+"/"+name+".fits"
        fits_with_copied_hdr(mask,main_image_path,data_object="Grid degraded to pixel",data_history="",fits_res_namepath=svpath,verbose=verbose)
        
    
    def get_subgrid_with_radius(self,r,r_type="pixel"):
        R = Length(r,setting=self.setting,input_type=r_type)
        if R.length_pix==self.Radius.length_pix:
            return self.get_circ_grid_flat()
        subgrid_select = [[gij.dist_from_center().length_pix<R.length_pix for gij in gi] for gi in self._sqrt_grid]
        subgrid_circ = self._flatten(self._sqrt_grid)[self._flatten(subgrid_select).tolist()]
        return subgrid_circ
    
    def get_subgrid_with_radius_from_nc(self,r,new_center,input_type="pixel"):
        R = Length(r,setting=self.setting,input_type=input_type)
        precision = self.precision_px if input_type=="pixel" else self.precision_radec
        C = Point(*new_center,[0,0],precision,self.setting,input_type)
        subgrid_select = [[ dist_btw_points(gij,C).length_pix<R.length_pix for gij in gi] for gi in self._sqrt_grid]
        subgrid_circ = self._flatten(self._sqrt_grid)[self._flatten(subgrid_select).tolist()]
        return subgrid_circ
    
    def _mask_center(self,grid_i,r_mask=0,r_type="pixel"):
        
        
        try:
            r_mask.length_pix
            R = r_mask
        except:
            R = Length(r_mask,setting=self.setting,input_type=r_type)
            
        
        if len(np.shape(grid_i))!=1:
            flat_grid_i = self._flatten(grid_i)
        else:
            flat_grid_i = grid_i
        grid_select    = [gi.dist_from_center().length_pix>R.length_pix for gi in flat_grid_i]
        masked_grid_i  = flat_grid_i[grid_select]
        return masked_grid_i