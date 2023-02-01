# create class combined_setting to be written by Multiplied_Df_post
# and read by H0_Combined(_reworked).py
# 17th Nov 2022

def test_zlist(zlist):
    if all(zlist!=z_list[0]):
        raise RuntimeError("The given redshifts are not identical. Something is off: z_list=",z_list)
    return 0

class combined_setting():
    def __init__(self,comment,z_source,z_lens,filters,setting_names):
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
        self.filters  = filters
        self.setting_names = setting_names
    
    def _test_zls(self):
        if self.z_lens>=self.z_source:
            raise RuntimeError("The redshift of the lens cannot be larger then or equal to the one of the source!")
        return 0