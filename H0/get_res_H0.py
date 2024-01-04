import os,sys
import pickle
import argparse as ap

from H0.tools_H0 import H0_Res
from Utils.statistical_tools import quantiles
from Utils.tools import present_program,success

def _get_H0_res(dir_name="PH0",res_dir="./results/",res_name="ph0_results.data"):
    PH0_resdir = res_dir+"/"+dir_name
    with open(f"{PH0_resdir}/{res_name}","rb") as f:
        PH0,H0 = pickle.load(f)
    return PH0,H0


def _print_H0_res(dir_name="PH0",res_dir="./results/",res_name="ph0_results.data"):
    PH0,H0  = _get_H0_res(dir_name=dir_name,res_dir=res_dir,res_name=res_name)
    h0_res,err_min,err_max = quantiles(PH0,H0,q=[0.16,.5,0.84],return_quantiles=False)
    H0_res = H0_Res(h0_res,[err_min,err_max])
    print("analytical: ", H0_res)

if __name__ == '__main__':
        
    present_program(sys.argv[0])
    
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Print res of H0",
                               formatter_class=ap.RawTextHelpFormatter)

    help_postdir   = "Name of the directory contains the H0 combined posterior"

    parser.add_argument("-dPH0","--dir_PostH0",dest='dir_ph0', type=str,default="PH0",
                        metavar='dir_ph0', action='store',
                        help=help_postdir)
    parser.add_argument("-rd","--res_dir",dest='res_dir', type=str,default="./results/",
                        metavar='res_dir', action='store',
                        help="Name of the results directory")
    parser.add_argument("-rn","--res_name",dest='res_name', type=str,default="ph0_results.data",
                        metavar='res_name', action='store',
                        help="Name of the results file")

    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")
    args      = parser.parse_args()
    dir_ph0   = args.dir_ph0
    res_dir   = args.res_dir
    res_name  = args.res_name
    verbose   = args.verbose
    _print_H0_res(dir_name=dir_ph0,res_dir=res_dir,res_name=res_name)
    success(sys.argv[0])
