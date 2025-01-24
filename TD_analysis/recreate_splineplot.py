import sys
import pycs3
import argparse as ap
import pathlib as pth

from Utils.tools import *
from TD_analysis.Standard_analysis_ABCD import set_orig_magshift, mask_A_peak


class ClassWrapper:
    def __init__(self, wrapped_class,str_new):
        self.wrapped_class = wrapped_class
        self._str = str_new

    def __str__(self):
        return self._str

    # Delegate all other attribute access and method calls to the wrapped class
    def __getattr__(self, attr):
        return getattr(self.wrapped_class, attr)

def recreate_splineplot(config,verbose=False):
    # copy from Standard_analysis to redo the analysis but only once
    config = get_config(config)
    lcs = config.get_lcs() 
    
    if config.maskA:
        lcs=mask_A_peak(lcs)
    
    
    savefig_path=pth.Path(config.analysis_directory)
    mkdir(savefig_path)
    
    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        if verbose:
            print("knt",knt)
        saveknt_path = savefig_path/str("analysis_kn"+str(knt))
        mkdir(savefig_path)        
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            if verbose:
                print("mltype",mltype_i)
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                mllist_name,mlfp = ml_config 
                if verbose:
                    print("ml_config",ml_config)
                # note: mlfp is for polyml the n* of Free Parameters ofthe polynomial, or the 
                kwargs_mc = {"knst_list":config.knotstep, "nit":1,
                 "dt_range":config.tsrand,"mc_res":config.mc_res}
                kwargs_ml = {"mltype":mltype_i,"mllist_name":mllist_name}

                if mltype_i=="polyml": 
                    kwargs_ml["mlfp"] = mlfp
                elif mltype_i=="splml": 
                    kwargs_ml["forcen"] = config.forcen
                    kwargs_ml["nmlspl"] = mlfp
                kwargs_mc.update(kwargs_ml)
                saveml_path      = saveknt_path/config.get_savemlpath(mltype_i,ml_config)# saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                savesplines_path = saveml_path/"splines/"
                
                try:
                    lcs, spline= pycs3.gen.util.readpickle( str(get_simdir(saveml_path))+'/initopt.pkl')
                    new_lcs = [] 
                    dt_corr = -2400000.5
                    for i,lc in enumerate(lcs):
                        lci = ClassWrapper(lc,str_new="[WST/"+ config.lcs_label[i]+"]")
                        lci.shifttime(dt_corr)
                        new_lcs.append(lci)

                    new_spline = ClassWrapper(spline,str_new="Intrinsic Spline")
                    new_spline.shifttime(dt_corr) 

                    nm=str(savesplines_path)+"/fit_lcs_w_spline.pdf"
                    print("Saving ",nm )
                    pycs3.gen.lc_func.display(new_lcs,[new_spline],nicefont=True,title="SDSSJ1433",showdates=True,filename=nm)
                except FileNotFoundError:
                    print("not modelled")

if __name__ == '__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Redo 1 analysis for each group (set of parameters for the lcs analysis) and redo the plot of the splines",
                               formatter_class=ap.RawTextHelpFormatter)
    help_lensname = "name of the lens to process"
    help_dataname = "name of the data set to process (Euler, SMARTS, ... )"
    parser.add_argument(dest='lensname', type=str,
                        metavar='lens_name', action='store',
                        help=help_lensname)
    parser.add_argument(dest='dataname', type=str,
                        metavar='dataname', action='store',
                        help=help_dataname)
    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")
    args = parser.parse_args()
    lensname = args.lensname
    dataname = args.dataname
    verbose  = args.verbose

    present_program(sys.argv[0])

    config = get_config(lensname=lensname,dataname=dataname,config_path="myconfig")
    recreate_splineplot(config)
    
    success(sys.argv[0])
