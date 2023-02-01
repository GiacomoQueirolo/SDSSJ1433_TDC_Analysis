import numpy as np
from copy import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

from tools import *
from conversion import conv_radec_to_xy,conv_xy_to_radec
from image_manipulation import get_rotangle,fits_with_copied_hdr, get_numPix


def dist(xy0,xy1):
    return np.sqrt( (xy0[0]-xy1[0])**2 + (xy0[1]-xy1[1])**2 ) 
    
def verify_type(_type):
    if _type not in ["pixel","radec"]:
        raise ValueError("Type either in pixel or radec, not "+str(_type))
        
class Point():
    def __init__(self,x,y,center,precision,setting,input_type="radec"):
        verify_type(input_type)
        self.setting = get_setting_module(setting,sett=True)
        if isinstance(center,Point):
            self.center_radec = center.get_radec()
            self.center_px    = center.get_pixel()
        else:        
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
        self.dist_from_center = self.get_dist_from_center()
        
    def get_kwargs(self,output_type="radec"):
        # for simplicity
        verify_type(output_type)
        if output_type=="radec":
            return {"center":self.center_radec,"precision":self.precision_radec,"setting":self.setting,"input_type":output_type}
        else:
            return {"center":self.center_px,   "precision":self.precision_px,   "setting":self.setting,"input_type":output_type}
    
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
        
    def get_dist_from_center(self):
        if not hasattr(self,"dist_from_center"):
            self.dist_from_center = Length(dist(self.center_radec,self.get_radec()),self.setting,input_type="radec")
        return  self.dist_from_center
    
    def dist_from(self,P2):
        return dist_btw_points(self,P2)
    
    def __add__(self,shift):
        kwarg_point = self.get_kwargs()
        if not isinstance(shift,Point) and type(shift,list):
            shift = Point(*shift,**kwarg_point)
            print("Warning: Assuming shift in radec")
        elif isinstance(shift,Point):
            if shift.center_radec!=self.center_radec:
                shift_center = [self.center_ra+shift.center_ra,self.center_dec+shift.center_dec]
                shift_ra  = shift.ra  + shift_center[0]
                shift_dec = shift.dec + shift_center[1]
                shift     = Point(shift_ra,shift_dec,**kwarg_point)
        shift = Shift(shift.ra,shift.dec,self.setting,"radec")
        ra  = self.ra+shift.dra
        dec = self.dec+shift.ddec
        return Point(ra,dec,**kwarg_point)

class Square(Point):
    def __init__(self,x,y,center,precision,setting,input_type="radec",radius=None):
        self.Point      = Point(x,y,center,precision,setting,input_type=input_type)
        super().__init__(x,y,**self.Point.get_kwargs(output_type=input_type))
        if not isinstance(center,Point):
            center = Point(*center,[0,0],precision,setting,input_type) 
        self.Center     = center
        self.area_px    = self.precision_px*self.precision_px
        self.area_radec = self.precision_radec*self.precision_radec
        if radius:
            if not isinstance(radius,Length):
                radius = Length(radius,precision,setting,input_type)
            self.inside_rad = self.dist_from_center<radius
                
    def __add__(self,shift):
        print("Warning: we are assuming this is a shift and moreover a shift of the whole grid")
        newPoint  = self.Point+shift
        newCenter = self.Center+shift
        newSquare = Square(newPoint.ra,newPoint.dec,newCenter,newPoint.precision_radec,
                        newPoint.setting,input_type="radec")
        # we assume is a shift of the whole grid: it doesn't change whether this is within or outside the radius
        if hasattr(self,"inside_rad"):
            newSquare.inside_rad = self.inside_rad
        return newSquare
    """
    def change_center(self,newCenter):
        if not isinstance(newCenter,Point):
            newCenter = Point(*newCenter,**self.get_kwargs())
            
        newSquare = copy(self)
        
        newSquare.center = newcenter
        return newSquare
    """

class Length():
    def __init__(self,length,setting,input_type="radec"): 
        self.setting = get_setting_module(setting,sett=True)
        if input_type == "radec":
            self.length_radec = length 
            self.length_pix   = length/self.setting.pix_scale
        elif input_type=="pixel":
            self.length_radec = length*self.setting.pix_scale
            self.length_pix   = length
    def __lt__(self,other): #less than
        return self.length_radec<other.lenght_radec
    def __le__(self,other): #less or equal
        return self.length_radec<=other.lenght_radec
    def __gt__(self,other): #greater than
        return self.length_radec>other.lenght_radec
    def __ge__(self,other): #less then
        return self.length_radec>=other.lenght_radec
    def __eq__(self,other): #equal
        return self.length_radec==other.lenght_radec
    def __ne__(self, other):#not equal
        return not(self.__eq__(self, other))
        
class Shift():
    def __init__(self,dx,dy,setting,input_type="radec"):
        if input_type=="radec":
            self.dra,self.ddec      = dx,dy
            self.dx_pix,self.dy_pix = np.reshape(conv_radec_to_xy(self.setting,self.dra,self.ddec),2).tolist()
        elif input_type=="pixel":
            self.dx_pix,self.dy_pix = dx,dy
            self.dra,self.ddec      = np.reshape(conv_xy_to_radec(self.setting,self.dx_pix,self.dy_pix),2).tolist()
        length = dist([dx,0],[0,dy])
        self.lenght = Length(length,setting,input_type)
        
def dist_btw_points(P1,P2):
    if not isinstance(P1,Point) or not isinstance(P2,Point):
        raise ValueError("The two imputs have to be instances of class Point. Else use dist(xy0,xy1)")
    if get_setting_name(P1.setting)!=get_setting_name(P2.setting):
        raise ValueError("The points given must have the same setting")
    ra1,dec1   = P1.get_radec()
    ra2,dec2   = P2.get_radec()
    dist_radec = dist([ra1,dec1],[ra2,dec2])
    return Length(dist_radec,P1.setting,input_type="radec")
    
class Square_Grid_Class():
    # create squared grid with centers of pixels
    # given center, "radius" (bc it can be made circular) 
    def __init__(self,center,radius,precision,setting,input_type="radec",r_mask=None,cnt_mask=None):
        self.Center = Point(*center,[0,0],precision,setting,input_type)
        self.Radius = Length(radius,setting,input_type)
        if r_mask:
            self.Rmask = Length(r_mask_cnt,setting,input_type)
            if not cnt_mask:
                self.Cmask = self.Center
            else:
                if isinstance(cnt_mask,Point):
                    self.Cmask = cnt_mask
                else:
                    self.Cmask = Point(*cnt_mask,[0,0],precision,setting,input_type) 
        else:
            self.Rmask = None
            self.Cmask = None
        self.precision_px    = self.Center.precision_px
        self.precision_radec = self.Center.precision_radec
        self.setting         = get_setting_module(setting,sett=True)
        self.points_radec    = self.initiate_radec_points()
        self.grid            = self.get_grid()
        if r_mask:
            self.mask()
        
    def _flatten(self,grid):
        try:
            return np.array(grid).reshape(len(grid)**2,2).tolist()
        except ValueError: 
            return np.array(grid).flatten().tolist()
            
    def initiate_radec_points(self):
        # let's first consider it in radec
        rd_radec = self.Radius.length_radec
        n  = int(2*rd_radec/self.precision_radec)
        n  = max([1,n]) 
        center_ra,center_dec = self.Center.ra,self.Center.dec
        
        dra  = np.linspace(center_ra-rd_radec,center_ra+rd_radec,n)
        ddec = np.linspace(center_dec-rd_radec,center_dec+rd_radec,n)
        points_centers_radec = [dra,ddec]
        points_radec         = self._flatten(np.transpose(np.meshgrid(*points_centers_radec)))
        return points_radec
        
    def get_grid(self):
        if hasattr(self,"grid"):
            return self.grid
        kwargs = {"center":self.Center,"precision":self.precision_radec,"setting":self.setting,\
                  "input_type":"radec","radius":self.Radius}
        grid = [Square(*gi,**kwargs) for gi in self.points_radec]
        return grid
        
    def _get_degraded_pixel_grid(self):
        Xc,Yc = self.Center.get_pixel()
        Rc    = self.Radius.length_pix
        topright_grid   = np.round([Xc+Rc,Yc+Rc],0).astype(int)
        bottomleft_grid = np.round([Xc-Rc,Yc-Rc],0).astype(int)
        
        xx = np.arange(bottomleft_grid[0],topright_grid[0])
        yy = np.arange(bottomleft_grid[1],topright_grid[1])
        
        points_pix = self.flatten(np.transpose(np.meshgrid(xx,yy)))
        
        if self.Rmask is not None:
            points_pix = self.mask(points_pix)
        return points_pix
        
    def _is_in_mask(self,point):
        if not isinstance(point,Point):
            point = Point(*point,self.center,self.precision,self.setting,self.input_type)
        return dist_btw_points(self.Cmask,point)<self.Rmask
        
        
    def mask(self,list_points=False):
        if list_points is False:
            self.grid = self.mask(self.grid)
            if hasattr(self,"degraded_pixel_grid"):
                self.degraded_pixel_grid = self.mask(self.degraded_pixel_grid)
        else:
            mask_points = [ not self._is_in_mask(list_points[i])  for i in range(len(mask_points)) ]
            list_points = np.array(list_points)[mask_points].tolist()
            return list_points
    
    def get_degraded_pixel_grid(self):
        if hasattr(self,"degraded_pixel_grid"):
            return self.degraded_pixel_grid
        self.degraded_pixel_grid = self._get_degraded_pixel_grid()
        return self.degraded_pixel_grid
    
    def get_area_square(self,output_type="radec"):
        verify_type(input_type)
        # they are the same for every square by construction
        if output_type=="pixel":
            return self.grid[0].area_px
        elif output_type=="radec":
            return self.grid[0].area_radec
            
    def get_totarea_grid(self,output_type="radec"):
        verify_type(input_type)
        # still the same area for each square
        return self.get_area_square(output_type)*len(self.get_grid())
    
    def plot_grid(self,color="green",ax=None,grid_to_plot="grid",alpha=.7):
        if not ax:
            fig,ax=plt.subplots()
        if grid_to_plot not in ["grid","degraded_pixel_grid"]:
            raise RuntimeError
        else:
            if grid_to_plot =="grid":
                class_point=True
                grid = self.get_grid()
                pix_width = self.precision_px
                angle     = get_rotangle(self.setting) 
                
            elif grid_to_plot=="degraded_pixel_grid":
                class_point=False # those are just coordinates
                grid = self.get_degraded_pixel_grid()
                pix_width = 1
                angle     = 0
        pix_height = pix_width
        PC = []
        for i in range(len(grid)):
            if not class_point:
                anchor_p = [grid[i][0]-pix_width/2.,grid[i][1]-pix_height/2.]
            else:
                anchor_p = [grid[i].x_pix-pix_width/2.,grid[i].y_pix-pix_height/2.]
            rect = patches.Rectangle(anchor_p, pix_width - 0.05*pix_width, pix_height - 0.05*pix_height,color = color,angle=angle)
            PC.append(rect)
        ax.add_collection(PatchCollection(PC,facecolor=color,alpha=alpha,edgecolor=None))
        # remember to do "ax.plot()" in order to plot it!
        return ax 
    
    def produce_degraded_fits(self,name="degraded_grid",verbose=True,_other_grid=None):
        # in order to be sure that the problem with the grid is not due to plotting
        # i now produce the fits file to compare w the original image
        # for the degraded pixel grid
        if not _other_grid:
            grid = self.get_degraded_pixel_grid_flat()
        else:
            grid = _other_grid
        main_image_path = self.setting.data_path+"/"+self.setting.image_name
        numpix=get_numPix(self.setting,twodim=True)
        mask = np.zeros(numpix)
        for px in grid:
            mask[px[0]][px[1]] = 1
        svpath = self.setting.data_path+"/"+name+".fits"
        fits_with_copied_hdr(mask,main_image_path,data_object="Grid degraded to pixel",data_history="",fits_res_namepath=svpath,verbose=verbose)
        
    
    def get_subgrid(self,r,input_type="pixel",newCenter=None):
        R = Length(r,setting=self.setting,input_type=r_type)
        if R==self.Radius:
            return self
        SubGrid = copy(self)
        if newCenter is None:
            newCenter = self.Center
        elif not isinstance(newCenter,"Point"):
            newCenter = Point(*newCenter,[0,0],self.precision,self.setting,input_type=input_type)
        Subgrid.Center = newCenter
        if R<self.Radius:
            Subgrid.Radius = R        
            subgrid_select = [ dist_btw_points(gi,C)<R for gi in self.grid]
            SubGrid.grid = np.array(self.grid)[subgrid_select].tolist()
            if hasattr(self,"degraded_pixel_grid"):
                subgrid_degraded_select = [ dist_btw_points(gi,C)<R for gi in self.degraded_pixel_grid]
                SubGrid.degraded_pixel_grid = np.array(self.degraded_pixel_grid)[subgrid_degraded_select].tolist()
            return SubGrid
        else:
            # pragma todo
            raise RuntimeError("Still to implement")
    
     
class Grid_Class(Square_Grid_Class):
    # create circular grid with centers of pixels
    # given center, radius
    def __init__(self,center,radius,precision,setting,input_type="radec",r_mask=None,cnt_mask=None):
        super.__init__(center,radius,precision,setting,input_type,r_mask,cnt_mask)
        self.orig_grid = copy(self.grid)
        self.circ_grid = self.get_circ_grid()
        if r_mask:
            self.mask()

    def get_circ_grid(self):
        if hasattr(self,"circ_grid"):
            return self.circ_grid
        grid_select = [gi.inside_rad for gi in self.orig_grid]
        grid_circ = np.array(self.orig_grid)[grid_select].tolist()
        if hasattr(self,"Rmask"):
            grid_circ = self.mask(grid_circ)
        return grid_circ
        
    def mask(self,list_points):
        if list_points is False:
            self.orig_grid = self.mask(self.orig_grid)
            self.circ_grid = self.mask(self.circ_grid)
            if hasattr(self,"degraded_pixel_grid"):
                self.degraded_pixel_grid = self.mask(self.degraded_pixel_grid)
        else:
            return super(Grid_Class, self).mask(list_points)
    
    def _get_degraded_pixel_circ_grid(self):
        sqrt_deg_pix = super(Grid_Class, self).get_degraded_pixel_grid()
        XYc = self.Center.get_pixel()
        Rc  = self.Radius
        grid_circ  = [ xyi for xyi in sqrt_deg_pix  if dist(xyi,XYc)<=Rc] 
        if self.Rmask is not None:
            grid_circ = self.mask(grid_circ)
        return grid_circ
    
    def get_degraded_pixel_circ_grid(self):
        if hasattr(self,"degraded_pixel_circ_grid"):
            return self.degraded_pixel_circ_grid
        self.degraded_pixel_circ_grid = self._get_degraded_pixel_circ_grid()
        return self.degraded_pixel_circ_grid
                
    def get_area_circ_grid(self,output_type="radec"):
        verify_type(output_type)
        # still the same area for each square
        return self.get_area_square(output_type)*len(self.get_circ_grid())
    
    def plot_grid_circ(self,color="green",ax=None,grid_to_plot="circ_grid",alpha=.7):
        if not ax:
            fig,ax=plt.subplots()
        if grid_to_plot not in ["circ_grid","degraded_pixel_circ_grid"]:
            raise RuntimeError
        else:
            if grid_to_plot =="circ_grid":
                class_point=True
                flat_grid = self.get_circ_grid()
                pix_width = self.precision_px
                angle     = get_rotangle(self.setting)            
            elif grid_to_plot=="degraded_pixel_circ_grid":
                grid        = self.get_degraded_pixel_circ_grid()
                class_point = False # those are just coordinates
                pix_width   = 1
                angle       = 0
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
    
    def produce_degraded_circ_fits(self,name="degraded_circ_grid",verbose=True):
        super(Grid_Class,self).produce_degraded_fits(name,verbose,_other_grid=self.get_degraded_pixel_circ_grid())
    
    def get_subgrid(self,r,input_type="pixel",newCenter=None):
        R = Length(r,setting=self.setting,input_type=r_type)
        if R==self.Radius:
            return self
        SubGrid = copy(self)
        if newCenter is None:
            newCenter = self.Center
        elif not isinstance(newCenter,"Point"):
            newCenter = Point(*newCenter,[0,0],self.precision,self.setting,input_type=input_type)
        Subgrid.Center = newCenter
        if R<self.Radius:
            Subgrid.Radius = R        
            subgrid_select = [ dist_btw_points(gi,C)<R for gi in self.circ_grid]
            SubGrid.circ_grid = np.array(self.circ_grid)[subgrid_select].tolist()
            if hasattr(self,"degraded_pixel_circ_grid"):
                subgrid_degraded_select = [ dist_btw_points(gi,C)<R for gi in self.degraded_pixel_circ_grid]
                SubGrid.degraded_pixel_circ_grid = np.array(self.degraded_pixel_circ_grid)[subgrid_degraded_select].tolist()
            return SubGrid
        else:
            # pragma TODO
            raise RuntimeError("Still to implement")
