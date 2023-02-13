from Utils.tools import *
from tabulate import tabulate

class Par():
    def __init__(self,sett,model,name, val,_min,_max,sig):
        self.sett = sett
        self.model = model
        self.name = name
        self.val = val
        self._min = _min
        self._max = _max
        self.sig  = sig
    def __str__(self):
        #return(self.name+": "+str(round(self.val,4))+"_-"+str(round(self._min,4))+"^+"+str(round(self._max,4))+"+-"+str(round(self.sig,4)))
        return(str(round(self.val,4))+"_-"+str(round(self._min,4))+"^+"+str(round(self._max,4))+"+-"+str(round(self.sig,4)))

if __name__=="__main__":
    k = ["settings_f105w_CPRP_ws.py",    "settings_f160w_CPRP.py",
         "settings_f140w_CPRP.py","settings_f475w_CPRP_ws.py",
         "settings_f160w_7030_CPRP.py","settings_f814w_CPRP_ws.py"]

    models     = ["Main_lens","Pert","Unif"]
    models_sub = models[1:]

    ignore_this_params = ["center_x","center_y","ra_image","dec_image","amp"]
    params = []
    for sett in k:
        s =get_setting_module(sett,1)
        #print("\n",strip_setting_name(s))
        llp = s.lens_light_params
        if s.sub:
            m = models_sub
        else:
            m=models
        for im,mod in enumerate(m):
            #print(mod)
            for prm in llp[0][im].keys():
                if prm in ignore_this_params:
                    continue
                #__print(params,llp[0][im][params],llp[-2][im][params],llp[-1][im][params], llp[1][im][params])
                par_i = Par(s,mod,prm,llp[0][im][prm],llp[-2][im][prm],llp[-1][im][prm], llp[1][im][prm])
                #print(par_i)
                params.append(par_i)
                
    
    #print(tabulate([['Alice', 24], ['Bob', 19]], headers=['Name', 'Age']))
    """
    Name      Age
    ------  -----
    Alice      24
    Bob        19"""

    #tab = [["",[Par_str for Par_sets in Par_sets ] ]]

    def empty_par(sett,model):
        return Par(sett,model,"",None,None,None,None)
        
    for mod  in models:
        print("\n")
        params_mod = [p if p.model==mod else empty_par(p.sett,mod) for p in params]
        tab = []
        names = list(set([pm.name for pm in params_mod]) )
        for nm in names:
            tab_i = [nm]
            values,mins,maxs,sigmas = [],[],[],[]
            for pm in params_mod:
                #if pm.name==nm:
                tab_i.append(str(pm))  
                values.append(pm.val)
                mins.append(pm._min)
                maxs.append(pm._max)
                sigmas.append(pm.sig)
            min_ = min(mins)
            max_ = max(maxs)
            val_ = .5*(max_+min_)
            sig_ = np.mean(sigmas)
            tab_i.append(str(Par("","","",val_,min_,max_,sig_)))
            tab.append(tab_i)
        print(tabulate(tab,headers=["Models Params",*[strip_setting_name(s) for s in k],"Combined"]))
