import sys
import pickle
import argparse
import numpy as np

from Utils.tools import *

def get_FR(Dmag):
    return (10**(-Dmag/2.5))
def get_sig_FR(dmag,sig_dmag):
    #return abs(np.log(10)*get_FR(dmag)*(-1/2.5)*sig_dmag)
    # wrong, bc it's non-linear and therefore the standard
    # error prop. function is not useful
    
    #print( 10**(-(dmag-sig_dmag)/2.5),get_FR(dmag))
    sig_max  = 10**(-(dmag-sig_dmag)/2.5)  -get_FR(dmag)
    sig_min  = get_FR(dmag) - 10**(-(dmag+sig_dmag)/2.5)
    sig  = (sig_max+sig_min)/2.

    return sig


#"J1433_forcen"
if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser        = argparse.ArgumentParser(description="Conversiong between Dmag to Flux Ratio with error propagation",formatter_class=argparse.RawTextHelpFormatter)
    help_lensname = "name of the lens to process"
    help_dataname = "name of the data set to process (Euler, SMARTS, ... )"

    parser.add_argument(dest='lensname', type=str,
                        metavar='lens_name', action='store',
                        help=help_lensname)
    parser.add_argument(dest='dataname', type=str,
                        metavar='dataname', action='store',
                        help=help_dataname)
    args = parser.parse_args()
    lensname = args.lensname
    dataname = args.dataname    
    config = lensname+"_"+dataname
    config = get_config(config)
    wddir = config.combined_directory 
    marginalisation_dir = wddir+"/"+ config.name_marg_spline +"/"
    comb_Dmag = pickle.load(open(marginalisation_dir  +"mag_"+ config.name_marg_spline + "_sigma_%2.2f" % config.sigmathresh + '_combined.pkl','rb'))
    from stnd_plot import print_res_w_err
    FR = [get_FR(res) for res in comb_Dmag.results]
    #print (comb_Dmag.results,comb_Dmag.tot_error)
    sig_FR = [get_sig_FR(res,err) for res,err in zip(comb_Dmag.results,comb_Dmag.tot_error)]
    for FRi,sig_FRi,nm in zip(FR,sig_FR,comb_Dmag.labels):
        print(FRi,"$\pm$",sig_FRi)
        print("FR ",nm,":",print_res_w_err(FRi,sig_FRi))
    success(sys.argv[0])
