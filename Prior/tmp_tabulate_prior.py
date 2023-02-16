import argparse
import numpy as np
from tabulate import tabulate

from Utils.tools import *
from Data.conversion import conv_radec_to_xy
from Data.input_data import init_lens_model_list
class Par():
    def __init__(self,sett,model,name, val,_min,_max,sig,latex=False):
        self.sett  = sett
        self.model = model
        self.name  = name
        self.val   = val
        self._min  = _min
        self._max  = _max
        self.sig   = sig
        self.latex = latex
    def __str__(self):
        """
        if not self.val is None:
            return(self.name+": "+str(round(self.val,4))+"_-"+str(round(self._min,4))+"^+"+str(round(self._max,4))+"+-"+str(round(self.sig,4)))
        else:
            return(self.name+": "+str(self.val)+"_-"+str( self._min )+"^+"+str( self._max)+"+-"+str(self.sig))
        """
        if self.latex:
            return(f'U({round(self._min,2) },{round(self._max,2)})')
        else:
            return(str(round(self.val,4))+" _"+str(round(self._min,4))+" ^"+str(round(self._max,4))+" +-"+str(round(self.sig,4)))

def wrap_Par(s,mod,i_mod,prm,kw_list,latex=False):
    return Par(sett=s,model=mod,name=prm,
                val=kw_list[0][i_mod][prm], 
                _min = kw_list[-2][i_mod][prm],
                _max=kw_list[-1][i_mod][prm],
                sig=kw_list[1][i_mod][prm],
                latex=latex)
                
def wrap_Par_ps(s,mod,i_mod,prm,i_ps,kw_list,latex=False):
    return Par(sett=s,model=mod,name=prm+"_"+str(i_ps),
                val=kw_list[0][i_mod][prm][i_ps], 
                _min = kw_list[-2][i_mod][prm][i_ps],
                _max=kw_list[-1][i_mod][prm][i_ps],
                sig=kw_list[1][i_mod][prm][i_ps],
                latex=latex)
                
def convert_par_radec_to_xy(par_ra,par_dec,approx=True,DS9=True):
    # give parameters in radec coord, convert it into pixel coords
    # -> problem: due to the rotation, the min-max are not 1-1 with the pixel coords
    # DS9 : add 1 to every coord, to be consistent with ds9
    ra,dec = par_ra.val,par_dec.val
    _min_ra,_min_dec = par_ra._min,par_dec._min
    _max_ra,_max_dec = par_ra._max,par_dec._max
    
    sig_ra,sig_dec = par_ra.sig,par_dec.sig
    x,y           =  conv_radec_to_xy(par_ra.sett,ra,dec)
    sig_x,sig_y   =  sig_ra/par_ra.sett.pix_scale,sig_dec/par_dec.sett.pix_scale 
    _edge_x11,_edge_y11 = conv_radec_to_xy(par_ra.sett,_min_ra,_min_dec)
    _edge_x12,_edge_y12 = conv_radec_to_xy(par_ra.sett,_min_ra,_max_dec)
    _edge_x21,_edge_y21 = conv_radec_to_xy(par_ra.sett,_max_ra,_min_dec)
    _edge_x22,_edge_y22 = conv_radec_to_xy(par_ra.sett,_max_ra,_max_dec)
    
    
    edge_x = [*_edge_x11,*_edge_x12,*_edge_x21,*_edge_x22]
    edge_y = [*_edge_y11,*_edge_y12,*_edge_y21,*_edge_y22]
    add_1 = 0
    if DS9:
        add_1 = 1
    x,y = x[0]+add_1,y[0]+add_1
    edge_x = np.array(edge_x)+add_1
    edge_y = np.array(edge_y)+add_1
    
    parx = Par(sett=par_ra.sett,model=par_ra.model,name=par_ra.name+"-> x",
                val=x,_min=None,_max =None,sig= sig_x)
    pary = Par(sett=par_dec.sett,model=par_dec.model,name=par_dec.name+"-> y",
                val=y,_min=None,_max =None,sig= sig_y)
    if approx:

        min_x = edge_x.min()
        max_x = edge_x.max()

        min_y = edge_y.min()
        max_y = edge_y.max()
        
        parx._max = max_x
        parx._min = min_x
        pary._max = max_y
        pary._min = min_y

        return parx,pary
    else:
        return parx,pary,edge_x,edge_y
        
def print_tabulate(settings,models,params,lens_light_models=False,only_diff=False):
    for mod  in models:
        print("\n")
        params_mod = [p for p in params  if p.model==mod]
        tab = []
        names = list(set([pm.name for pm in params_mod]) )
        for nm in names:
            tab_i = [nm]
            values,mins,maxs,sigmas = [],[],[],[]
            for pm in params_mod:
                if pm.name==nm:
                    tab_i.append(str(pm))  
                    values.append(pm.val)
                    mins.append(pm._min)
                    maxs.append(pm._max)
                    sigmas.append(pm.sig)
            if only_diff:
                # print only the differences   
                if all(mins[0]==mn for mn in mins) and all(maxs[0]==mx for mx in maxs) and all(values[0]==val for val in values) and all(sigmas[0]==sg for sg in sigmas):
                    continue

            min_ = min(mins)
            max_ = max(maxs)
            val_ = .5*(max_+min_)
            sig_ = np.mean(sigmas)
            tab_i.append(str(Par("","","",val_,min_,max_,sig_)))
            tab.append(tab_i)
        s_k = settings
        if lens_light_models:
            if mod==models[0]:
                s_k = [ strip_setting_name(s) for s in settings if not get_setting_module(s,1).sub]
        print(tabulate(tab,headers=[mod,*s_k,"Combined"]))



def _better_prm_name(name):
    if "_" in name:
        if "center_" in name:
            return name[-1].capitalize()
        nm_spl = name.replace("_","_{")
        nmspl = name.split("_")
        if len(nmspl)==2:
            name = "".join([nmspl[0],"_{",nmspl[1].capitalize(),"}"])
        else:
            pass
    if "e1"==name or "e2"==name:
        return "e_"+str(name[-1])
    if "psi" in name or "gamma" in name:
        return "\\"+name
    return name
    
def print_tab_latex(settings,models,params,lens_light_models=False):
    for mod  in models:
        print("\n")
        params_mod = [p for p in params  if p.model==mod]
        tab = []
        names = list(set([pm.name for pm in params_mod]) )
        for nm in names:
            tab_i = [_better_prm_name(nm)]
            for pm in params_mod:
                if pm.name==nm:
                    tab_i.append(str(pm))
            tab.append(tab_i)
            
        s_k = [get_filter(s) for s in settings]
        if lens_light_models:
            if mod==models[0]:
                s_k = [ get_filter(s) for s in settings if not get_setting_module(s,1).sub]
        print(tabulate(tab,headers=[mod,*s_k],tablefmt="latex_raw"))




if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Print table of prior value (tmp)")
    parser.add_argument('-od','--only_diff',action="store_true", dest="only_diff", default=False,help="Only print different parameters")
    parser.add_argument('-l','--latex',action="store_true", dest="latex", default=False,help="Print the LateX table")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args          = parser.parse_args()
    settings      = args.SETTING_FILES
    only_diff     = args.only_diff
    latex         = args.latex
    models_ll     = ["Main_lens","Pert","Unif"]
    models_ll_sub = models_ll[1:]
    model_lm      = [init_lens_model_list(s) for s in settings]
    model_source  = ["Sersic"]
    model_ps      = ["Images pos"]
    
    if all(mlm==model_lm[0] for mlm in model_lm):
        model_lm = ["Main Lens ("+model_lm[0][0]+")","Pert ("+model_lm[0][1]+")","Shear ("+model_lm[0][2]+")"]
    else:
        raise RuntimeError("Comparing different mass models!") 

    ignore_this_params = ["amp"]#["center_x","center_y","ra_image","dec_image","amp"]
    
    params_ll  = [] # lens light
    params_lm  = [] # lens mass
    params_src = [] # source light
    params_psx = [] # point sources x
    params_psy = [] # point sources y
    for sett in settings:
        s =get_setting_module(sett,1)
        #print("\n",strip_setting_name(s))
        llp = s.lens_light_params
        llm = s.lens_params
        ps  = s.ps_params
        if not s.WS:
            sp = s.source_params
        if s.sub:
            mll = models_ll_sub
        else:
            mll = models_ll
            
        for im,mod in enumerate(mll):
            for prm in llp[0][im].keys():
                if prm in ignore_this_params:
                    continue
                par_i = wrap_Par(s,mod,im,prm,llp,latex=latex)
                params_ll.append(par_i)
        
        for im, mod in enumerate(model_lm):
            for prm in llm[0][im].keys():
                if prm in ignore_this_params:
                    continue
                par_i = wrap_Par(s,mod,im,prm,llm,latex=latex)
                params_lm.append(par_i)

        for im,mod in enumerate(model_source):
            if s.WS:
                continue
            for prm in sp[0][im].keys():
                if prm in ignore_this_params:
                    continue
                par_i = wrap_Par(s,mod,im,prm,sp,latex=latex)
                params_src.append(par_i)
                
        for i_mod,mod in enumerate(model_ps):

            for prm in ps[0][i_mod].keys():
                if prm in ignore_this_params:
                    continue
                for i in range(len(ps[0][i_mod][prm])):
                    par_i = wrap_Par_ps(s,mod,i_mod,prm,i,ps,latex=latex) 
                    if prm=="ra_image":
                        params_psx.append(par_i)
                    elif prm=="dec_image":
                        params_psy.append(par_i)
                
        

    #print(tabulate([['Alice', 24], ['Bob', 19]], headers=['Name', 'Age']))
    """
    Name      Age
    ------  -----
    Alice      24
    Bob        19"""    
    # point sources pos
    print("\nPoint Sources params")
    print("####################")
    if latex:
        print_tab_latex(settings,model_ps,params_psx)
        print_tab_latex(settings,model_ps,params_psy)
    else:
        print_tabulate(settings,model_ps,params_psx,lens_light_models=False,only_diff=only_diff)
        print_tabulate(settings,model_ps,params_psy,lens_light_models=False,only_diff=only_diff)
    
    # lens mass params
    print("\nLens mass params")
    print("################")
    if latex:
        print_tab_latex(settings,model_lm,params_lm)
    else:
        print_tabulate(settings,model_lm,params_lm,lens_light_models=False,only_diff=only_diff)
       
    # lens light params
    print("\nLens light params")
    print("#################")
    if latex:
        print_tab_latex(settings,models_ll,params_ll,lens_light_models=True)
    else:
        print_tabulate(settings,models_ll,params_ll,lens_light_models=True,only_diff=only_diff)
        
    # source light params
    if params_src!=[]:
        print("\nSource light params")
        print("###################")
        settings = list(set([get_setting_name(p.sett) for p in params_src]))
        if latex:
            print_tab_latex(settings,model_source,params_src)
        else:
            print_tabulate(settings,model_source,params_src,lens_light_models=False,only_diff=only_diff)
