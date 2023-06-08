import numpy as np
import matplotlib.pyplot as plt

from Utils.tools import get_setting_module,get_savefigpath,check_if_SUB
from Data.input_data import init_kwrg_data
from Data.image_manipulation import get_lens_light_model
from Plots.plotting_tools import plot_model_WS as PM
from Plots.plotting_tools import get_ModelPlot
from Plots.PSO_plotlight import get_vlim
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Data.conversion import conv_xy_to_radec,conv_radec_to_xy

def coordinate_arrows(setting,ax, frame_size, color='w', font_size=15, arrow_size=0.05):
    d  = frame_size
    d0 = d / 8.
    p0 = d / 15.
    pt = d / 9.
    deltaPix = setting.pix_scale # coords.pixel_width
    ra0, dec0 = np.reshape(conv_xy_to_radec(setting,(d - d0) / deltaPix, d0 / deltaPix),2)#coords.map_pix2coord((d - d0) / deltaPix, d0 / deltaPix)
    xx_, yy_ =  np.reshape(conv_radec_to_xy(setting,ra0,dec0),2) #coords.map_coord2pix(ra0, dec0)
    xx_ra, yy_ra = np.reshape(conv_radec_to_xy(setting,ra0 + p0, dec0),2)#coords.map_coord2pix(ra0 + p0, dec0)
    xx_dec, yy_dec = np.reshape(conv_radec_to_xy(setting,ra0, dec0 + p0),2)#coords.map_coord2pix(ra0, dec0 + p0)
    xx_ra_t, yy_ra_t = np.reshape(conv_radec_to_xy(setting,ra0+pt,dec0),2)#coords.map_coord2pix(ra0 + pt, dec0)
    xx_dec_t, yy_dec_t = np.reshape(conv_radec_to_xy(setting,ra0,dec0+pt),2)#coords.map_coord2pix(ra0, dec0 + pt)

    ax.arrow(xx_ * deltaPix, yy_ * deltaPix, (xx_ra - xx_) * deltaPix, (yy_ra - yy_) * deltaPix,
             head_width=arrow_size * d, head_length=arrow_size * d, fc=color, ec=color, linewidth=1)
    ax.text(xx_ra_t * deltaPix, yy_ra_t * deltaPix, "E", color=color, fontsize=font_size, ha='center')
    ax.arrow(xx_ * deltaPix, yy_ * deltaPix, (xx_dec - xx_) * deltaPix, (yy_dec - yy_) * deltaPix,
             head_width=arrow_size * d, head_length=arrow_size * d, fc
             =color, ec=color, linewidth=1)
    ax.text(xx_dec_t * deltaPix, yy_dec_t * deltaPix, "N", color=color, fontsize=font_size, ha='center')


print("Reconstructing the 'results_model' plot, but considering the lens light model when this is previously subtracted")

#sett="f140w_CPRP_PLLI"
sett="f814w_CPRP_PLLI_ws"
setting         = get_setting_module(sett,1)
modelPlot       = get_ModelPlot(setting) 
savefig_path    = get_savefigpath(sett)
v_min,v_max     = setting.v_min,setting.v_max
res_min,res_max = setting.res_min,setting.res_max

#PM(modelPlot,savefig_path,v_min,v_max,res_min,res_max,band=(0,"tmp"))
kwdata,mask = init_kwrg_data(setting,return_mask=True)
Orig        =  mask*kwdata["image_data"]
modelPlot_i = modelPlot._band_plot_list[0]
CL_Conv     = modelPlot_i._model
if check_if_SUB(setting):
    # add the lens light model   
    ll_model = get_lens_light_model(setting)
    CL_Conv = CL_Conv + ll_model
    if setting.already_sub:
        Orig += mask*ll_model 
kw_vlim  = get_vlim(Orig,fact=7)
Res      = modelPlot_i._norm_residuals
frame_size= max(np.shape(Res))*5
extnt = [0, frame_size, 0, frame_size]
fontsize=15
fig,ax = plt.subplots(1,3,figsize=(16, 8))
imgs = [Orig,CL_Conv,Res]
titles = ["Original","Reconstructed","Norm. Residual"]

for i in range(3):
    c_vertical = 1/13. #+ fontsize / d / 10.**2
    c_horizontal = 1./30
    color='w'

    if i!=2:
        im = ax[i].imshow(imgs[i],origin="lower", cmap = plt.cm.gist_heat,extent=extnt,**kw_vlim )#vmin=setting.v_min,vmax=setting.v_max)
        colorbar_label=r'log$_{10}$ Flux'
        ax[i].text(frame_size * c_horizontal, frame_size - frame_size * c_vertical, titles[i],color="w", backgroundcolor="k",fontsize=fontsize)
    else:
        im=ax[i].imshow(imgs[i],origin="lower", cmap='bwr',extent=extnt,vmin=setting.res_min,vmax=setting.res_max)
        colorbar_label=r'(f${}_{\rm model}$ - f${}_{\rm data}$)/$\sigma$'
        ax[i].text(frame_size * c_horizontal, frame_size - frame_size * c_vertical, titles[i],color="k", backgroundcolor="w",fontsize=fontsize)
        color="k"

    divider = make_axes_locatable(ax[i])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(im, cax=cax)    
    cb.set_label(colorbar_label,
                 fontsize=fontsize)

    ax[i].get_xaxis().set_visible(False)
    ax[i].get_yaxis().set_visible(False)
    ax[i].autoscale(False)
    coordinate_arrows(setting,ax[i], frame_size, color=color,arrow_size=0.02, font_size=fontsize)
    p0   = frame_size / 15.
    dist = (frame_size/max(np.shape(Orig)))*1/setting.pix_scale
    ax[i].plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
    ax[i].text(p0 + dist / 2., p0 + 0.01 * frame_size, '1"', fontsize=fontsize, color=color, ha='center')
savename = f"{get_savefigpath(setting)}/reconstructed_results_model.pdf"
plt.tight_layout()
plt.savefig(savename)
print(f"Saved {savename}")

