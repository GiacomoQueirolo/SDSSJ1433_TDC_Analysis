# presubtract the modelled lensed source light for F160W

from tools import *
from create_fits import source_light 
from image_manipulation  import load_fits,fits_with_copied_hdr
import sys
import numpy as np
from argparse import ArgumentParser



if __name__=="__main__":
    present_program(sys.argv[0])
    parser = ArgumentParser(description="Create fits file with all light components")
    parser.add_argument('SETTING_FILE_REF',default=None,help="Reference setting file for F160W to consider")
    parser.add_argument('SETTING_FILE_OUT',default=None,help="Ouput setting file for F160W to consider")
    args = parser.parse_args()
    setting_ref = get_setting_module(args.SETTING_FILE_REF,sett=True)
    setting_out = get_setting_module(args.SETTING_FILE_OUT,sett=True)
    SL_Conv = source_light(setting_ref,unconvolved=False)
    
    in_namepath = setting_ref.data_path+"/"+setting_ref.image_name
    Image   = load_fits(in_namepath)
    
    Image_SUBS = Image - SL_Conv
    out_namepath = setting_out.data_path+"/"+setting_out.image_name
    history = str(in_namepath+" - lensed source light model from "+strip_setting_name(setting_ref))
    fits_with_copied_hdr(Image_SUBS,in_namepath,data_history=history,fits_res_namepath=out_namepath,overwrite=True,verbose=True) 
    success(sys.argv[0])