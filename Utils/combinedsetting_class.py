# create class combined_setting to be written by Multiplied_Df_post
# and read by H0_Combined(_reworked).py
# 17th Nov 2022

from os import getcwd 
from tools import get_filter,get_combined_setting_module

def test_zlist(zlist):
    if any([z!=zlist[0] for z in zlist]):
        raise RuntimeError("The given redshifts are not identical. Something is off: z_list=",zlist)
    return 0

class combined_setting():
    def __init__(self,comment,z_source,z_lens,filters,setting_names,savedir,KDE=False,BC=True):
        self.comment  = comment
        if type(z_source) is not list:
            self.z_source = z_source
        else:
            test_zlist(z_source)
            self.z_source = z_source[0]
        if type(z_lens) is not list:
            self.z_lens   = z_lens
        else:
            test_zlist(z_lens)
            self.z_lens = z_lens[0]
        self._test_zls()
        self.filters = filters
        self.savedir = savedir
        self.KDE     = KDE
        self.BC      = BC   # if True, the Df combinations are AB,AC,BC, else AB,AC,AD
        if self.BC:
            from Posterior_analysis.fermat_pot_analysis import labels_Df_BC as labels_Df
            from Posterior_analysis.mag_remastered import labels_Rmag_BC as labels_Rmag
        else:
            from Posterior_analysis.fermat_pot_analysis import labels_Df_AD as labels_Df
            from Posterior_analysis.mag_remastered import labels_Rmag_AD as labels_Rmag
        self.labels_Df     = labels_Df      
        self.labels_Rmag   = labels_Rmag     
        self.setting_names = setting_names
        # absolute path for the main lenstronomy directory
        # this only works if the program that calls this class was run in
        # the right directory. This should be always true, else we better hard-code the path 
        self.main_dir_abs = getcwd()


    def _test_zls(self):
        if self.z_lens>=self.z_source:
            raise RuntimeError("The redshift of the lens cannot be larger then or equal to the one of the source!")
        return 0
    
    def get_respath(self):
        # return absolute path to the result directory
        return f"{self.main_dir_abs}/{self.savedir}/"
        

    def gen_cmb_sett_name(self,cmb_sett_name=None):
        if hasattr(self,"cmb_sett_name"):
            self.cmb_sett_name = self.check_cmb_sett_name(self.cmb_sett_name)
            if cmb_sett_name is not None:
                if cmb_sett_name!=self.cmb_sett_name:
                    raise RuntimeError(f"This combined setting already have a name: {self.cmb_sett_name}, I can't change it to {cmb_sett_name}")
            return self.cmb_sett_name
        if cmb_sett_name is None:
            dir_name      = self.savedir.split("/")[-2]
            cmb_sett_name = "cmb_setting_"+[get_filter(s)+"_" for s in self.setting_names]+dir_name
        cmb_sett_name = self.check_cmb_sett_name(cmb_sett_name)
        print("Combined_setting_name: ",cmb_sett_name)
        self.cmb_sett_name = cmb_sett_name
        return self.cmb_sett_name

    def check_cmb_sett_name(self,cmb_sett_name):
        if not "cmb_setting" in cmb_sett_name:
            cmb_sett_name = "cmb_setting_"+cmb_sett_name 
        try:
            #_ = find_combined_setting_path(cmb_sett_name)
            #raise FileExistsError(f"This name for combined_setting_name \"{cmb_sett_name}\" already exists. Try another")
            # this is not necessarily a problem
            # it is if it doen't have the same settings
            old_cmb_sett = get_combined_setting_module(cmb_sett_name)
            if self!=old_cmb_sett:
                raise FileExistsError(f"This name for combined_setting_name \"{cmb_sett_name}\" already exists but have different settings. Find another name")
            else:
                print(f"This name for combined_setting_name \"{cmb_sett_name}\" already exists but it's the same. Overwriting it")
                return cmb_sett_name
        except FileNotFoundError:
            pass
        return cmb_sett_name
    
    def __eq__(self,other_cmb_sett):
        if set(self.setting_names)==set(other_cmb_sett.setting_names) and set(self.filters)==set(other_cmb_sett.filters):
            if self.z_lens==other_cmb_sett.z_lens and self.z_source==other_cmb_sett.z_source:
                if self.KDE==other_cmb_sett.KDE:
                    return True
        return False
