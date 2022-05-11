# define a conversion of the name of each paramater for the plots
import re

def standard(param):
        UOM="[\"]"
        param_title = param

        if param=="theta_E_lens0":
            param_title=r"$\theta_E^1$"
            UOM="[\"]"
        elif param=="theta_E_lens1":
            param_title=r"$\theta_E^2$"
            UOM="[\"]"
        elif param=="e1_lens0":
            param_title=r"e$_1$"
            UOM="[]"
        elif param=="e2_lens0":
            param_title=r"e$_2$" 
            UOM="[]"
        elif param=="center_x_lens0":
            param_title=r"x$_{lens}$"
            UOM="[\"]"
        elif param=="center_y_lens0":
            param_title=r"y$_{lens}$"
            UOM="[\"]"
        elif param=="R_sersic_lens_light0":
            param_title=r"R$_{Sersic}^{pert. light}$"
        elif param=="n_sersic_lens_light0":
            param_title=r"n$_{Sersic}^{pert. light}$"
            UOM="[]"
        elif param=="R_sersic_source_light0":
            param_title=r"R$_{Sersic}^{source light}$"
        elif param=="n_sersic_source_light0":
            param_title=r"n$_{Sersic}^{source light}$"
            UOM="[]"
        elif param=="e1_source_light0":
            param_title=r"e$_{1}^{source light}$"
            UOM="[]"
        elif param=="e2_source_light0" :
            param_title=r"e$_{2}^{source light}$"
            UOM="[]"
        elif param=="center_x_lens_light0":
            param_title=r"x$_{pert.}$"
            UOM="[\"]"
        elif param=="center_y_lens_light0":
            param_title=r"y$_{pert.}$"
            UOM="[\"]"
        elif "gamma_ext" in param:
            param_title=r"$\gamma_{Shear}$"
            UOM=""
        elif "psi_ext" in param:
            param_title=r"$\psi_{Shear}$"
            UOM="[°]"    
        else:
            raise TypeError("Not found "+param)
        return param_title,UOM
        
def with_main_lens_light(param):
    if not re.search(".+_lens_light.+",param): 
        return standard(param)
    else:
        UOM="[\"]"
        if param=="center_x_lens_light0":
            param_title=r"x$_{lens}$"        
        elif param=="center_y_lens_light0":
            param_title=r"y$_{lens}$"
        elif param=="R_sersic_lens_light0":
            param_title=r"R$_{Sersic}^{lens light}$"
        elif param=="n_sersic_lens_light0":
            param_title=r"n$_{Sersic}^{lens light}$"
            UOM="[]"
        elif param=="e1_lens_light0":
            param_title=r"$e_1^{lens light}$"
            UOM="[]"
        elif param=="e2_lens_light0":
            param_title=r"$e_2^{lens light}$"
            UOM="[]"
        elif param=="center_x_lens_light1":
            param_title=r"x$_{pert.}$"
        elif param=="center_y_lens_light1":
            param_title=r"y$_{pert.}$"
        elif param=="R_sersic_lens_light1":
            param_title=r"R$_{Sersic}^{pert. light}$"
        elif param=="n_sersic_lens_light1":
            param_title=r"n$_{Sersic}^{pert. light}$"
            UOM="[]"
        else:
            raise TypeError("Not found "+param)
        return param_title,UOM     
    
def newProfile(param):
    if not "gamma_lens" in param:
        return standard(param)
    else:
        param_title=r"$\gamma_1$"
        UOM = "[\"]"
        return param_title, UOM

def newProfile_noPert(param):
    UOM = "[\"]"
    if param=="center_x_lens1" :
        param_title=r"x$_{pert.}$"
        return param_title, UOM    
    elif param=="center_y_lens1":
        param_title=r"y$_{pert.}$"
        return param_title, UOM    
    else:
        return newProfile(param)
        
def newProfile_with_main_lens_light(param):
    """UOM = "[\"]"
    if param=="center_x_lens_light0" :
        param_title=r"x$_{lens}$"
        return param_title, UOM    
    elif param=="center_y_lens_light0" :
        param_title=r"y$_{lens}$"
        return param_title, UOM 
    elif param=="center_x_lens_light1" :
        param_title=r"x$_{pert.}$"
        return param_title, UOM    
    elif param=="center_y_lens_light1" :
        param_title=r"y$_{pert.}$"
        return param_title, UOM 
    else:
        return newProfile(param)"""
    if not re.search(".+_lens_light.+",param): 
        return newProfile(param)
    else:
        UOM="[\"]"
        if param=="center_x_lens_light0":
            param_title=r"x$_{lens}$"        
        elif param=="center_y_lens_light0":
            param_title=r"y$_{lens}$"
        elif param=="R_sersic_lens_light0":
            param_title=r"R$_{Sersic}^{lens light}$"
        elif param=="n_sersic_lens_light0":
            param_title=r"n$_{Sersic}^{lens light}$"
            UOM="[]"
        elif param=="e1_lens_light0":
            param_title=r"$e_1^{lens light}$"
            UOM="[]"
        elif param=="e2_lens_light0":
            param_title=r"$e_2^{lens light}$"
            UOM="[]"
        elif param=="center_x_lens_light1":
            param_title=r"x$_{pert.}$"
        elif param=="center_y_lens_light1":
            param_title=r"y$_{pert.}$"
        elif param=="R_sersic_lens_light1":
            param_title=r"R$_{Sersic}^{pert. light}$"
        elif param=="n_sersic_lens_light1":
            param_title=r"n$_{Sersic}^{pert. light}$"
            UOM="[]"
        else:
            raise TypeError("Not found "+param)
        return param_title,UOM     
