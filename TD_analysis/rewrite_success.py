import sys
import argparse as ap
import pathlib as pth

from Utils.tools import *
from Standard_analysis_ABCD import check_analysis


if __name__ == '__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Rewrite/update success file",
                               formatter_class=ap.RawTextHelpFormatter)
    help_lensname = "name of the lens to process"
    help_dataname = "name of the data set to process"
    parser.add_argument(dest='lensname', type=str,
                        metavar='lens_name', action='store',
                        help=help_lensname)
    parser.add_argument(dest='dataname', type=str,
                        metavar='dataname', action='store',
                        help=help_dataname)
    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")
    args      = parser.parse_args()
    lensname  = args.lensname
    dataname  = args.dataname
    verbose   = args.verbose
    name_prog = sys.argv[0]
    present_program(name_prog)

    sys.path.append("myconfig/")
    config_file = "myconfig_" + lensname + "_" + dataname
    config = importlib.import_module(config_file)
    config = get_config(config)
    savefig_path=pth.Path(config.analysis_directory)
    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("analysis_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                mllist_name,mlfp = ml_config 
                saveml_path      = saveknt_path/config.get_savemlpath(mltype_i,ml_config)# saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                with open(str(saveml_path)+"/td.data", 'rb') as f:
                    timedelays = pickle.load( f)
                kw_success = check_analysis(timedelays,max_std=5,second_chance=True)
                with open(str(saveml_path)+"/success.data","wb") as f:
                    pickle.dump(kw_success,f)
    success(sys.argv[0])