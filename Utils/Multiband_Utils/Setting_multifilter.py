from Utils.Multiband_Utils.tools_multifilter import _append
from Utils.Multiband_Utils.tools_multifilter import *

from Utils.tools import *
from Utils.get_res import load_whatever

# util program to define the setting class

# standard multiflter directory
stnd_mutlifilter_setting_dir = "multifilter_setting"
class Setting_Multifilter():
    def __init__(self,settings,backup_path,multifilter_setting_dir=stnd_mutlifilter_setting_dir,name=None,main_dir="./"):
        self.backup_path = backup_path
        self.settings    = get_setting_module(settings,1)
        self.main_dir    = main_dir
        self.setting_dir = create_path_from_list([main_dir,multifilter_setting_dir])

        if name is None:
            name      = "_".join([strip_setting_name(s) for s in self.settings])
        self.name     = name
        self.savename = verify_multifilter_setting_name(self.name)
        self._create_savefigdir()
        self._create_savemcmcdir()
        self.backend_filename = create_path_from_list([self.savemcmc_path,f"backend_mcmc_{name}.hp5"])
        self._check_z()
        self._check_if_CP()
        self._check_if_fixedmag()
        self.prior_lens = _check_same_prior_lens(self.settings,ret_prior=True)
        self.pll        = self.prior_lens.pll
        self.filters  = [s.filter_name for s in self.settings]
        self.comments = create_comments(self.settings)
        self.allWS    = all([check_if_WS(sett)  for sett in self.settings ])
        self.allSUB   = all([check_if_SUB(sett) for sett in self.settings ])
        self.arrange_params_lens()
        self.arrange_params_lens_light()
        self.arrange_params_ps()
        if not self.allWS:
            self.arrange_params_source()
        self._savemyself()
        pass
    def _create_savefigdir(self):
        self.savefig_path = create_path_from_list([self.backup_path,self.name])
        mkdir(self.savefig_path)
    def _create_savemcmcdir(self):
        self.savemcmc_path = create_path_from_list([self.backup_path,"mcmc_"+self.name])
        mkdir(self.savemcmc_path)
    def _check_z(self):
        self.z_lens   =  self.settings[0].z_lens
        self.z_source =  self.settings[0].z_source
        for sett in self.settings[1:]:
            if sett.z_lens   != self.z_lens:
                raise RuntimeWarning("The given setting files have different z_lens")
            if sett.z_source != self.z_source:
                raise RuntimeWarning("The given setting files have different z_lens")
        return 0
    def _check_if_CP(self):
        self.CP = check_if_CP(self.settings[0])
        for sett in self.settings:
            if check_if_CP(sett) !=self.CP:
                raise RuntimeError("The given settings have different model (CP)")
        return 0
    def _check_if_fixedmag(self):
        self.fixed_mag = self.settings[0].fixed_mag
        for sett in self.settings:
            if sett.fixed_mag !=self.fixed_mag:
                raise RuntimeError("The given settings have choice overt the fixed magnification")
        return 0

    def arrange_params_lens(self):
        # the lens parameters should all be the same apart from the coordinates
        # we assume the frame of reference is correct enough
        # and give the first as parameter for the lens
        self.lens_params = self.settings[0].lens_params

    def arrange_params_ps(self):
        kwargs_ps_init = []
        kwargs_ps_sigma = []
        fixed_ps = []
        kwargs_lower_ps = []
        kwargs_upper_ps = []
        for sett in self.settings:
            _append(kwargs_ps_init, sett.ps_params[0])
            _append(kwargs_ps_sigma, sett.ps_params[1])
            _append(fixed_ps, sett.ps_params[2])
            _append(kwargs_lower_ps, sett.ps_params[3])
            _append(kwargs_upper_ps, sett.ps_params[4])
        self.ps_params = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]

    def arrange_params_lens_light(self):
        kwargs_lens_light_init = []
        kwargs_lens_light_sigma = []
        fixed_lens_light = []
        kwargs_lower_lens_light = []
        kwargs_upper_lens_light = []
        
        for sett in self.settings:
            """
            # this should be the lenght of the modelled profiles
            if check_if_SUB(sett):                
                _append(kwargs_lens_light_init, [])
                _append(kwargs_lens_light_sigma, [])
                _append(fixed_lens_light, [])
                _append(kwargs_lower_lens_light, [])
                _append(kwargs_upper_lens_light, [])
            """

            _append(kwargs_lens_light_init, sett.lens_light_params[0])
            _append(kwargs_lens_light_sigma, sett.lens_light_params[1])
            _append(fixed_lens_light, sett.lens_light_params[2])
            _append(kwargs_lower_lens_light, sett.lens_light_params[3])
            _append(kwargs_upper_lens_light, sett.lens_light_params[4])
        
        self.lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]

    def arrange_params_source(self):
        kwargs_source_init = []
        kwargs_source_sigma = []
        fixed_source = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        
        for sett in self.settings:
            # this should be the lenght of the modelled profiles
            if check_if_WS(sett):
                continue
                """
                _append(kwargs_source_init, [])
                _append(kwargs_source_sigma, [])
                _append(fixed_source,[])
                _append(kwargs_lower_source,[])
                _append(kwargs_upper_source, [])
                """
            _append(kwargs_source_init, sett.source_params[0])
            _append(kwargs_source_sigma, sett.source_params[1])
            _append(fixed_source, sett.source_params[2])
            _append(kwargs_lower_source, sett.source_params[3])
            _append(kwargs_upper_source, sett.source_params[4])
            
                
        self.source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    
    def get_savejson_name(self,filename):
        savejson_name = self.savename.replace("mltf_setting",filename).replace(".dll",".json")
        return savejson_name
    def get_savejson_path(self,filename): 
        savejson_path = self.savemcmc_path+"/"+ self.get_savejson_name(filename=filename)
        return savejson_path
    
    def savejson_data(self,data,filename):
        savepath= self.get_savejson_path(filename)
        save_json(data=data,filename=savepath)

    def get_data(self,filename):
        # note: to be general, the file has already to be defined with the correct end
        filepath= self.savemcmc_path+"/"+filename
        return load_whatever(filepath)

    def get_mcmc(self):
        kw_mcmc = {}
        for kw in ["smpl","prm","logL"]:
            mcmc_kw = "mcmc_"+kw
            name_mcmc = self.get_savejson_name(mcmc_kw)
            if kw=="prm":
                name_mcmc.replace("json","dat")
            kw_mcmc[mcmc_kw] = self.get_data(name_mcmc)
        return kw_mcmc

    def _savemyself(self):
        with open(f"{self.setting_dir}/{self.savename}","wb") as f:
            dill.dump(self,f)
            
    def __str__(self) -> str:
        return self.name

@check_setting
def _check_same_prior_lens(settings,ret_prior=False):
    prior_lens = settings[0].lens_prior
    if len(settings)==1:
        if ret_prior:
            return prior_lens
        return 0
    for sett in settings[1:]:
        if prior_lens!=sett.lens_prior:
            raise RuntimeWarning("Lens prior is not consistent for all setting files")
    if ret_prior:
        return prior_lens
    return 0

@check_setting
def create_multifilter_setting(settings,backup_path,multifilter_setting_dir=None,name=None,main_dir="./"):
    settings = get_setting_module(settings,1)
    if multifilter_setting_dir is None:
        print(f"Defining multifilter setting directory as the standard one: {stnd_mutlifilter_setting_dir}")
        multifilter_setting_dir = stnd_mutlifilter_setting_dir
    return Setting_Multifilter(settings=settings,backup_path=backup_path,name=name,\
                               multifilter_setting_dir=multifilter_setting_dir,main_dir=main_dir)

def create_comments(settings):
    cmnt= "Multifilter setting file considering setting files:\n"+" ".join([get_setting_name(sett)for sett in settings ])
    cmnt+="\nwith the following comments:"
    for sett in settings:
        cmnt+="\n"+sett.comments
    return cmnt
