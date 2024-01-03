# define a conversion of the name of each paramater for the plots
import re

def standard(param):
        UOM="[\"]"
        param_title = param

        if param=="theta_E_lens0":
            param_title=r"$\theta_\mathrm{E}^\mathrm{ML}$"
            UOM="[\"]"
        elif param=="theta_E_lens1":
            param_title=r"$\theta_\mathrm{E}^{\mathrm{Pert}}$"
            UOM="[\"]"
        elif param=="e1_lens0":
            param_title=r"e$_1^\mathrm{ML}$"
            UOM="[]"
        elif param=="e2_lens0":
            param_title=r"e$_2^\mathrm{ML}$" 
            UOM="[]"
        elif param=="center_x_lens0":
            param_title=r"x$^\mathrm{ML}$"
            UOM="[\"]"
        elif param=="center_y_lens0":
            param_title=r"y$^\mathrm{ML}$"
            UOM="[\"]"
        elif param=="R_sersic_lens_light0":
            param_title=r"R$_{Sersic}^{Pert light}$"
        elif param=="n_sersic_lens_light0":
            param_title=r"n$_{Sersic}^{Pert light}$"
            UOM="[]"
        elif param=="R_sersic_source_light0":
            param_title=r"R$_{Sersic}^{Host light}$"
        elif param=="n_sersic_source_light0":
            param_title=r"n$_{Sersic}^{Host light}$"
            UOM="[]"
        elif param=="e1_source_light0":
            param_title=r"e$_{1}^{Host light}$"
            UOM="[]"
        elif param=="e2_source_light0" :
            param_title=r"e$_{2}^{Host light}$"
            UOM="[]"
        elif param=="center_x_lens_light0":
            param_title=r"x$^{\mathrm{Pert}}$"
            UOM="[\"]"
        elif param=="center_y_lens_light0":
            param_title=r"y$^{\mathrm{Pert}}$"
            UOM="[\"]"
        elif "gamma_ext" in param:
            param_title=r"$\gamma^\mathrm{Shear}$"
            UOM="[]"
        elif "psi_ext" in param:
            param_title=r"$\psi^\mathrm{Shear}$"
            UOM=r"[$^\circ$]"    
        else:
            raise TypeError("Not found "+param)
        return param_title,UOM
        
def with_main_lens_light(param):
    if not re.search(".+_lens_light.+",param): 
        return standard(param)
    else:
        UOM="[\"]"
        if param=="center_x_lens_light0":
            param_title=r"x$^{\mathrm{ML}}$"        
        elif param=="center_y_lens_light0":
            param_title=r"y$^{\mathrm{ML}}$"
        elif param=="R_sersic_lens_light0":
            param_title=r"R$_{\mathrm{Sersic}}^{\mathrm{ML light}}$"
        elif param=="n_sersic_lens_light0":
            param_title=r"n$_{\mathrm{Sersic}}^{\mathrm{ML light}}$"
            UOM="[]"
        elif param=="e1_lens_light0":
            param_title=r"$e_1^{\mathrm{ML light}}$"
            UOM="[]"
        elif param=="e2_lens_light0":
            param_title=r"$e_2^{\mathrm{ML light}}$"
            UOM="[]"
        elif param=="center_x_lens_light1":
            param_title=r"x$^{\mathrm{Pert}}$"
        elif param=="center_y_lens_light1":
            param_title=r"y$^{\mathrm{Pert}}$"
        elif param=="R_sersic_lens_light1":
            param_title=r"R$_{\mathrm{Sersic}}^{\mathrm{Pert light}}$"
        elif param=="n_sersic_lens_light1":
            param_title=r"n$_{\mathrm{Sersic}}^{\mathrm{Pert light}}$"
            UOM="[]"
        else:
            raise TypeError("Not found "+param)
        return param_title,UOM     
   

def newProfile(param):
    if not "gamma_lens" in param:
        return standard(param)
    else:
        param_title=r"$\gamma^\mathrm{ML}$"
        UOM = "[\"]"
        return param_title, UOM

def newProfile_noPert(param):
    UOM = "[\"]"
    if param=="center_x_lens1" :
        param_title=r"x$^{Pert}$"
        return param_title, UOM    
    elif param=="center_y_lens1":
        param_title=r"y$^{Pert}$"
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
        param_title=r"x$^{Pert}$"
        return param_title, UOM    
    elif param=="center_y_lens_light1" :
        param_title=r"y$^{Pert}$"
        return param_title, UOM 
    else:
        return newProfile(param)"""
    if not re.search(".+_lens_light.+",param): 
        return newProfile(param)
    else:
        UOM="[\"]"
        if param=="center_x_lens_light0":
            param_title=r"x$^\mathrm{ML}$"        
        elif param=="center_y_lens_light0":
            param_title=r"y$^\mathrm{ML}$"
        elif param=="R_sersic_lens_light0":
            param_title=r"R$_{Sersic}^{ML light}$"
        elif param=="n_sersic_lens_light0":
            param_title=r"n$_{Sersic}^{ML light}$"
            UOM="[]"
        elif param=="e1_lens_light0":
            param_title=r"$e_1^{ML light}$"
            UOM="[]"
        elif param=="e2_lens_light0":
            param_title=r"$e_2^{ML light}$"
            UOM="[]"
        elif param=="center_x_lens_light1":
            param_title=r"x$^{Pert}$"
        elif param=="center_y_lens_light1":
            param_title=r"y$^{Pert}$"
        elif param=="R_sersic_lens_light1":
            param_title=r"R$_{Sersic}^{Pert light}$"
        elif param=="n_sersic_lens_light1":
            param_title=r"n$_{Sersic}^{Pert light}$"
            UOM="[]"
        else:
            raise TypeError("Not found "+param)
        return param_title,UOM
        
def newProfile_FS(param):
    if not "center_" in param and not "source_light" in param:
        return newProfile(param)
    elif param=="center_x_source_light0":
        param_title=r"$x^{Host}$"
        UOM = "[\"]"
    elif param=="center_y_source_light0":
        param_title=r"$y^{Host}$"
        UOM = "[\"]"
    else:
        raise TypeError("Not found "+param)
    return param_title, UOM

def newProfile_w_images(param):
    if not "image" in param:
        return newProfile(param)
    else:
        kw_im = {"0":"A","1":"B","2":"C","3":"D"}
        UOM = "[\"]"
        radec = param.split("_")[0]
        nradec = param[-1]
        param_title = r"${"+radec+r"}_{"+kw_im[nradec]+r"}$"
        return param_title, UOM
