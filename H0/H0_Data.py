# Tools to streamline the work when calling Df and Dt posteriors
import os,dill
from datetime import datetime
from Utils.get_res import get_combined_Df
from TD_analysis.pycs_get_res import get_combined_res
from Utils.tools import mkdir,get_config,get_combined_setting_respath,get_config_respath,save_log_command


def create_PH0resdir(res_dir="./results/",dir_ph0="PH0/",verbose=False,overwrite=False):
    mkdir(res_dir)
    PH0_resdir = res_dir+"/"+dir_ph0
    if os.path.isfile(PH0_resdir):
        if os.listdir(PH0_resdir):
            if verbose:
                print(f"Previous results found in {PH0_resdir}")
            if overwrite:
                if verbose:
                    print(f"Overwriting them (actually moving them to old_{PH0_resdir})")
                os.rename(PH0_resdir, "old_"+PH0_resdir)
            else:
                td = datetime.today().strftime("%y%m%d")
                if verbose:
                    print(f"Changing result directory name: {PH0_resdir}->{PH0_resdir}_{td}")
                PH0_resdir = PH0_resdir+"_"+td
    print("Saving results in: ",PH0_resdir)
    mkdir(PH0_resdir)
    return PH0_resdir

def create_links(Df_or_Dt,dir_src,dir_dest,overwrite=False):
    # dir_src and dir_dest should be the absolute path
    try:
        os.symlink(dir_src,dir_dest)
    except FileExistsError:
        if overwrite:
            print(f"Warning: Symlink found for {Df_or_Dt} posterior. This should not happen, but I am instructed to overwrite it")
            os.unlink(dir_dest)
            os.symlink(dir_src,dir_dest)
        else:
            raise RuntimeError(f"Previous symlink found for  {Df_or_Dt} posterior, this should not happen")
    return 0

def get_Dt_post(dt_name,PH0_resdir="results/PH0",pycs_path="./my_pycs_scripts/",configpath="/myconfig/",link=True,overwrite=False):
    dt_conf = get_config(dt_name,config_path=configpath,main_dir=pycs_path)
    Dt_res  = get_combined_res(dt_conf,main_dir_path=pycs_path)
    
    dir_dt    = dt_conf.get_respath()
    dt_resdir = PH0_resdir+"/Dt_post"
    if link:
        create_links("Dt",dir_src=dir_dt,dir_dest=dt_resdir,overwrite=overwrite)
        # saving config file
        with open(f"{PH0_resdir}/dt_config.dll","wb") as f:
            dill.dump(dt_conf,f)
    print("For time delay using posterior : ",str(dt_name))
    return Dt_res

def get_Df_post(df_name,PH0_resdir="results/PH0",lenstronomy_path="./lenstronomy/",link=True,overwrite=False):
    combined_setting, Combined_PDF,Combined_bins_or_points = get_combined_Df(df_name,main_dir=lenstronomy_path)
    dir_df    = combined_setting.get_respath()
    df_resdir = PH0_resdir+"/Df_post"
    if link:
        create_links("Df",dir_src=dir_df,dir_dest=df_resdir,overwrite=overwrite)
        # saving combined setting
        with open(f"{PH0_resdir}/df_cmb_sett.dll","wb") as f:
            dill.dump(combined_setting,f)
    print("For fermat pot. using posterior: ",str(dir_df))
    return combined_setting,Combined_PDF,Combined_bins_or_points
    
def write_readme(path,dt_name,df_name,lenstronomy_path="./lenstronomy/",pycs_path="./my_pycs_scripts/",configpath="/myconfig/",log_command=True):
    # save a txt file in the final path specifing the path
    # used for to find the dt and df data for the posterior

    df_path          = get_combined_setting_respath(df_name,main_dir=lenstronomy_path)
    dt_path          = get_config_respath(dt_name,config_path=configpath,main_dir=pycs_path)
    with open(str(path)+"/paths_to_data_used.txt","w") as f:
        f.write("Time delay posterior used:\n")
        f.write(str(dt_path)+"\n")
        f.write("Fermat potential posterior used:\n")
        f.write(str(df_path)+"\n")
    if log_command:
        save_log_command(path)
    return 0
