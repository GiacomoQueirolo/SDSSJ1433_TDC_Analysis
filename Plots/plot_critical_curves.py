#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Utils.tools import *
from Utils.get_res import *
from Data.input_data import *
from Data.image_manipulation import *

import sys
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

from lenstronomy.Data.pixel_grid import PixelGrid
from lenstronomy.LensModel.lens_model_extensions import LensModelExtensions


# In[ ]:


def plot_line_set(ax, coords, line_set_list_x, line_set_list_y, origin=None, flipped_x=False,  
                  pixel_offset=True, *args, **kwargs):
    """
    plotting a line set on a matplotlib instance where the coordinates are defined in pixel units with the lower left
    corner (defined as origin) is by default (0, 0). The coordinates are moved by 0.5 pixels to be placed in the center
    of the pixel in accordance with the matplotlib.matshow() routine.

    :param ax: matplotlib.axis instance
    :param coords: Coordinates() class instance
    :param origin: [x0, y0], lower left pixel coordinate in the frame of the pixels
    :param line_set_list_x: numpy arrays corresponding of different disconnected regions of the line
     (e.g. caustic or critical curve)
    :param line_set_list_y: numpy arrays corresponding of different disconnected regions of the line
     (e.g. caustic or critical curve)
    :param color: string with matplotlib color
    :param flipped_x: bool, if True, flips x-axis
    :param points_only: bool, if True, sets plotting keywords to plot single points without connecting lines
    :param pixel_offset: boolean; if True (default plotting), the coordinates are shifted a half a pixel to match with
     the matshow() command to center the coordinates in the pixel center
    :return: plot with line sets on matplotlib axis in pixel coordinates
    """
    if origin is None:
        origin = [0, 0]
    pixel_width = coords.pixel_width
    pixel_width_x = pixel_width
    if 'linestyle' not in kwargs:
        kwargs['linestyle'] = ""
    if 'marker' not in kwargs:
        kwargs['marker'] = "o"
    if 'markersize' not in kwargs:
        kwargs['markersize'] = 1
    if flipped_x:
        pixel_width_x = -pixel_width
    if pixel_offset is True:
        shift = 0.5
    else:
        shift = 0
    for i in range(len(line_set_list_x)):
        x_c, y_c = coords.map_coord2pix(line_set_list_x[i], line_set_list_y[i])
        ax.plot((x_c + shift) * pixel_width_x + origin[0], (y_c + shift) * pixel_width + origin[1],*args, **kwargs)
    return ax

def critical_curves_plot(ax, pixel_grid, lens_model, kwargs_lens,  color_crit='r',
                    pixel_offset=False,return_radec=False, *args, **kwargs):
    """
    :param ax: matplotlib axis instance
    :param pixel_grid: lenstronomy PixelGrid() instance (or class with inheritance of PixelGrid()
    :param lens_model: LensModel() class instance
    :param kwargs_lens: lens model keyword argument list
    :param color_crit: string, color of critical curve
    :param pixel_offset: boolean; if True (default plotting), the coordinates are shifted a half a pixel to match with
     the matshow() command to center the coordinates in the pixel center
    :param args: argument for plotting curve
    :param kwargs: keyword arguments for plotting curves
    :return: updated matplotlib axis instance
    """

    lens_model_ext = LensModelExtensions(lens_model)
    #pixel_width = pixel_grid.pixel_width
    #frame_size = np.max(pixel_grid.width)
    #coord_center_ra, coord_center_dec = pixel_grid.center
    kw_pg = kw_pixel_grid(pixel_grid)
    #ra_crit_list, dec_crit_list = lens_model_ext.critical_curve_tiling(kwargs_lens, compute_window=frame_size,
    #                                                                     start_scale=pixel_width, max_order=10,
    #                                                                     center_x=coord_center_ra,
    #                                                                     center_y=coord_center_dec)
    ra_crit_list, dec_crit_list = lens_model_ext.critical_curve_tiling(kwargs_lens,  max_order=10, **kw_pg)
    ra0, dec0 = pixel_grid.radec_at_xy_0
    origin = [ra0, dec0]
    if return_radec:
        kwrg_crit_curve = {"ra":ra_crit_list,"dec":dec_crit_list,"pixel_grid":pixel_grid,"origin":origin} 
        return kwrg_crit_curve


    plot_line_set(ax, pixel_grid, ra_crit_list, dec_crit_list, color=color_crit, origin=origin,
                             pixel_offset=False, *args,
                            **kwargs)
    return ax
    
def kw_pixel_grid(pixel_grid):
    start_scale = pixel_grid.pixel_width
    compute_window = np.max(pixel_grid.width)
    center_x, center_y = pixel_grid.center
    #ra0, dec0 = pixel_grid.radec_at_xy_0
    #origin = [ra0, dec0]
    kw = {"start_scale":start_scale,"compute_window":compute_window,"center_x":center_x,"center_y":center_y}
    return kw
        

def _caustic_plot(ra_points,dec_points,kwargs_lens,lens_model):
    ra_caustics, dec_caustics = lens_model.ray_shooting(ra_points, dec_points, kwargs_lens)
    return ra_caustics, dec_caustics

def caustic_plot(ax, pixel_grid, lens_model, kwargs_lens,  color_crit='r',
                    pixel_offset=False,return_radec=False, *args, **kwargs):
    """
    :param ax: matplotlib axis instance
    :param pixel_grid: lenstronomy PixelGrid() instance (or class with inheritance of PixelGrid()
    :param lens_model: LensModel() class instance
    :param kwargs_lens: lens model keyword argument list
    :param color_crit: string, color of critical curve
    :param pixel_offset: boolean; if True (default plotting), the coordinates are shifted a half a pixel to match with
     the matshow() command to center the coordinates in the pixel center
    :param args: argument for plotting curve
    :param kwargs: keyword arguments for plotting curves
    :return: updated matplotlib axis instance
    """

    lens_model_ext = LensModelExtensions(lens_model)
    pixel_width = pixel_grid.pixel_width
    frame_size = np.max(pixel_grid.width)
    coord_center_ra, coord_center_dec = pixel_grid.center
    ra0, dec0 = pixel_grid.radec_at_xy_0
    origin = [ra0, dec0]
    ra_crit_list, dec_crit_list = lens_model_ext.critical_curve_tiling(kwargs_lens, compute_window=frame_size,
                                                                         start_scale=pixel_width, max_order=10,
                                                                         center_x=coord_center_ra,
                                                                         center_y=coord_center_dec)
    
    ra_caustic_list,dec_caustic_list = _caustic_plot(ra_crit_list,dec_crit_list,kwargs_lens,lens_model)
    if return_radec:
        return ra_caustic_list ,dec_caustic_list

    plot_line_set(ax, pixel_grid, ra_caustic_list,dec_caustic_list, color=color_crit, origin=origin,
                             pixel_offset=False, *args,**kwargs)
    return ax

    
def transf_fits(sett,fits,multi=1):
    new_trsnfmtrx = get_transf_matrix(sett,in_arcsec=False)/multi
    fits.header["CD1_1"]=new_trsnfmtrx[0][0]
    fits.header["CD1_2"]=new_trsnfmtrx[0][1]
    fits.header["CD2_1"]=new_trsnfmtrx[1][0]
    fits.header["CD2_2"]=new_trsnfmtrx[1][1]
    coord_A    = [218.3449834 , 60.12145745] # ra,dec
    coord_A_xy = newpixel_grid.map_coord2pix(sett.ps_params[0][0]["ra_image"][0],sett.ps_params[0][0]["dec_image"][0])
    fits = shift_astrometry(fits,coord_A,coord_A_xy)
    return fits

def plot_cross(fits_data,xy_pix,edge_len=10,multi=1):
    edge_len = int(edge_len*multi)
    x,y = int(xy_pix[0]),int(xy_pix[1])
    for i in range(edge_len):
        fits_data[x][y+i] = 1
        fits_data[x][y-i] = 1
        fits_data[x+i][y] = 1
        fits_data[x-i][y] = 1
    return fits_data


# In[ ]:


if __name__=="__main__":
    present_program(sys.argv[0])
    parser = ArgumentParser(description="Plot critical curves")
    parser.add_argument("-sf","--singlefile",dest="singlefile", default=False,action="store_true",
                        help="Produce only 1 FITS file")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    
    args       = parser.parse_args()
    settings   = get_setting_module(args.SETTING_FILES,sett=True)
    multi      = 5
    fits       = True
    singlefile = args.singlefile
    for sett in settings:
        print(get_setting_name(sett))
        kwres_lens = get_kwres(sett)["kwargs_results"]["kwargs_lens"]
        lens_model = init_lens_model(sett)
        nx,ny      = get_numPix(sett,twodim=True)
        trsnf_mtrx = get_transf_matrix(sett,in_arcsec=True)
        pixel_grid = PixelGrid(nx,ny,trsnf_mtrx,sett.ra_at_xy_0,sett.dec_at_xy_0)
        
        fig,ax     = plt.subplots(1,1,figsize=(12,12))     
        #ax_cc      = critical_curves_plot(ax, pixel_grid, lens_model, kwres_lens,  color_crit='r', pixel_offset=False,return_radec=False)#,*args, **kwargs)
        kw_cc      = critical_curves_plot(ax, pixel_grid, lens_model, kwres_lens,  color_crit='r', pixel_offset=False,return_radec=True)
        ra_crit,dec_crit = kw_cc["ra"],kw_cc["dec"]
        ra_caus,dec_caus = _caustic_plot(ra_points=ra_crit,dec_points=dec_crit,kwargs_lens=kwres_lens,lens_model=lens_model)
        
        kwres_images = get_kwres(sett)["kwargs_results"]["kwargs_ps"][0]
        radec_images = [[ra,dec] for ra,dec in zip(kwres_images["ra_image"],kwres_images["dec_image"])]
        
        kwres_source = load_whatever(get_savefigpath(sett)+"/read_source.data")
        radec_source = kwres_source["source_ra"][0],kwres_source["source_dec"][0]
        
        if fits:
            canvas = np.zeros((nx*multi,ny*multi))
            if not singlefile:
                canvas2 = np.zeros((nx*multi,ny*multi))
            newpixel_grid = PixelGrid(nx*multi,ny*multi,trsnf_mtrx/multi,sett.ra_at_xy_0,sett.dec_at_xy_0)
            for i in range(len(ra_crit)):
                x_cc, y_cc = newpixel_grid.map_coord2pix(ra_crit[i],dec_crit[i])
                x_cc, y_cc = int(x_cc),int(y_cc)
                canvas[y_cc][x_cc] = 1
            for i in range(len(ra_caus)):
                x_cs, y_cs = newpixel_grid.map_coord2pix(ra_caus[i],dec_caus[i])
                x_cs, y_cs = int(x_cs),int(y_cs)
                if singlefile:
                    canvas[y_cs][x_cs] = 1
                else:
                    canvas2[y_cs][x_cs] = 1
                    
            # add cross at images pos and source pos 
            images = [newpixel_grid.map_coord2pix(radec_i[0],radec_i[1]) for radec_i in radec_images]
            source = newpixel_grid.map_coord2pix(*radec_source)
            for i in range(len(images)):
                canvas = plot_cross(canvas,images[i])
            canvas2    = plot_cross(canvas2,source)

            image_file = sett.data_path+sett.image_name
            fits_path  = get_savefigpath(sett)+"/crit_curves_"+sett.filter_name+".fits"
            
            fits       = fits_with_copied_hdr(data=canvas,data_object="Critical Curves",fits_parent_path=image_file,fits_res_namepath=None)
            fits       = transf_fits(sett,fits,multi)
            fits.writeto(fits_path, overwrite=True)
            print("Saving "+fits_path)
            if not singlefile:
                fits_path2 = get_savefigpath(sett)+"/caustics_"+sett.filter_name+".fits"
                fits2 = fits_with_copied_hdr(data=canvas2,data_object="Caustics",fits_parent_path=image_file,fits_res_namepath=None)
                fits2 = transf_fits(sett,fits2,multi)            
                fits2.writeto(fits_path2, overwrite=True)
                print("Saving "+fits_path2)   
        else:
            origin = list(pixel_grid.radec_at_xy_0)
            ax_CC     = plot_line_set(ax, pixel_grid, comb_ra, comb_dec, origin=origin)
            fig.savefig(get_savefigpath(sett)+"/crit_curves.png")
            print("Saving "+get_savefigpath(sett)+"/crit_curves.png")
    success(sys.argv[0])

    

