# Copy from HST_HR_Main.py - we test the MOD_LLFR
# from 26th oktober 22 
# add a new logL term to consider the mag ratio


import os,sys
import corner
import importlib
import numpy as np
import json,copy,pickle
from astropy.io import fits
from datetime import datetime
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from lenstronomy.Data.psf import PSF
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.PointSource.point_source import PointSource


# In[2]:


import get_res
from tools import *
import mag_remastered
from image_manipulation import *
# MOD_LLFR
from custom_logL import logL_combined
from pycs_get_res import get_combined_mag
from my_pycs_scripts.from_Dmag_to_FR import get_FR,get_sig_FR

from input_data import init_lens_model,init_kwrg_data,init_kwrg_psf,init_kwrg_numerics,get_kwargs_model


# In[ ]:


############################
present_program(sys.argv[0])
############################


# In[ ]:


parser = ArgumentParser(description="Lens modelling program")
parser.add_argument('-rt','--run_type',type=int,dest="run_type",default=0,help="Type of run: \n \
                0 = standard, PSO_it = 400*rf   PSO_prt = 200*rf      MCMCr = 800*rf  MCMCb = 200*rf\n \
                1 = append  MCMCr=800*rf\n\
                2 = test run  PSO_it = 3   PSO_prt = 3   MCMCr = 2 MCMCb = 1\n\
                3 = append test   MCMCr=1\n\
                (PSO_it: PSO iterations, PSO_prt: PSO particles, MCMCr: MCMC run steps, MCMCb: MCMC burn in steps)")
parser.add_argument('-rf','--run_factor',type=float,dest="run_factor",default=20.,help="Run factor to have longer run")
parser.add_argument('-tc','--threadCount',type=int,dest="threadCount",default=150,help="Number of CPU threads for the MCMC parallelisation (max=160)")
parser.add_argument('SETTING_FILE',default="",help="setting file to model")

args         = parser.parse_args()
run_type     = args.run_type
run_fact     = args.run_factor
setting_name = get_setting_name(args.SETTING_FILE).replace(".py","")
threadCount  = args.threadCount  
RND = False #set a random start of the PSO
n_run_cut = 50  # to re-implement
#Model PSO/MCMC settings
append_MC=False
if run_type==0:
    n_iterations = int(400*run_fact) #number of iteration of the PSO run
    n_particles  = int(200*run_fact) #number of particles in PSO run
    n_run  = int(800*run_fact) #MCMC total steps 
    n_burn = int(200*run_fact) #MCMC burn in steps
elif run_type ==1:
    append_MC=True
    n_run  = int(800*run_fact) #MCMC total steps 
elif run_type==2:
    n_iterations = int(3) #number of iteration of the PSO run
    n_particles  = int(3) #number of particles in PSO run
    n_run  = int(2) #MCMC total steps 
    n_burn = int(1) #MCMC burn in steps
elif run_type==3:
    append_MC   = True
    n_run  = int(1) #MCMC total steps 
else:
    raise RuntimeError("Give a valid run_type or implement it your own")

np.seterr(all="ignore");


# In[3]:


backup_path   = "backup_results"
savemcmc_path = get_savemcmcpath(setting_name,backup_path)
savefig_path  = get_savefigpath(setting_name,backup_path) 
mkdir(savefig_path)
mkdir(savemcmc_path)
setting_path = find_setting_path(setting_name)
os.system("cp "+setting_path+"/"+setting_name+".py "+savefig_path+".") #we copy the setting file to that directory
sys.path.append(setting_path)


# In[6]:


setting = get_setting_module(setting_name).setting()
CP = check_if_CP(setting)
WS = check_if_WS(setting)


# In[ ]:


if CP:
    print("WARNING: Considering the PEMD mass profile for the main lens")
if WS:
    print("WARNING: this model DO NOT consider the source")


# In[ ]:


if WS and setting_name[-3:]!="_ws":
    setting_dir=find_setting_path(setting_name)
    new_name = setting_dir+"/"+setting_name+"_ws.py"
    os.system("cp "+setting_dir+"/"+setting_name+".py "+new_name )
    def line_prepender(filename, line):
        with open(filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\r\n') + '\n' + content)
    line_prepender(new_name,"#SIMPLE COPY FROM "+setting_name+".py")
    setting_name=setting_name+"_ws"


# In[ ]:


# datetime object containing current date and time
now = datetime.now()
# dd/mm/YY H:M:S
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

print("Running: ",sys.argv[0])
print("Machine: ",os.uname()[1]) 
print("Setting file:", setting_name)
print("Started the :", dt_string)


# In[ ]:


"""
image_file   = setting.data_path+setting.image_name
err_file     = setting.data_path+setting.err_name
psf_file     = setting.data_path+setting.psf_name 
err_psf_file = setting.data_path+setting.err_psf
"""

##Printout of the results
print_res = open(savefig_path+"results.txt","w")
print_res.write("Results for "+setting_name+" \n")
print_res.write("#################################\n\n")
print_res.write("Date and time of start : "+ dt_string+"\n")
print_res.write("Results obtained with "+sys.argv[0]+"\n")
print_res.write("Setting file used settings/"+setting_name+"\n")
print_res.write("Comments: "+setting.comments+"\n")
print_res.write("append_MC: "+str(append_MC)+"\n")


# In[ ]:


"""image = load_fits(image_file)
plot_image(image,setting,savefig_path+"/original_image.png")
numPix = image.shape[0]
err_image = load_fits(err_file)
# the masking procedure define the mask, subtract the eventual lens light profile,
# plot the error and image masked and return the subtracted (unmasked) image and the mask
mask = create_mask(image,setting)
# Subtract lens light (if necessary)
image_sub = subtract_light(image,setting)
# Plot the subtracted and masked data image and error
if setting.sub:
    plot_image(image_sub*mask,setting,savefig_path+"/image_sub_mask.png")
else:
    plot_image(image_sub*mask,setting,savefig_path+"/image_mask.png")
plot_image(err_image*mask,setting,savefig_path+"/err_mask.png",err_image=True)"""


# In[ ]:


"""kwargs_data = { 'noise_map':   err_image,  # noise map
                'ra_at_xy_0':  setting.ra_at_xy_0,  # RA  at (0,0) pixel
                'dec_at_xy_0': setting.dec_at_xy_0, # DEC at (0,0) pixel 
                'ra_shift' :   0., #shifts the coordinate system with respect to 'ra_at_xy_0'    
                'dec_shift':   0., #shifts the coordinate system with respect to 'dec_at_xy_0'    
                'transform_pix2angle': setting.transform_pix2angle,#Mpix2coord,  
                # matrix to translate shift in pixel in shift in relative RA/DEC (2x2 matrix). 
                # Make sure it's units are arcseconds or the angular units you want to model.
                'image_data': image_sub  # 2d data vector
              }"""
kwargs_data,mask = init_kwrg_data(setting,saveplots=True,backup_path=backup_path,return_mask=True)


# ### PSF

# In[ ]:


"""psf_image = load_fits(psf_file)

if setting.pssf>1:
    plot_image(psf_image,setting,savefig_path+"/psf_supersampled.png")
else:
    plot_image(psf_image,setting,savefig_path+"/psf.png")

#We import the psf error image 
err_psf_image = load_fits(err_psf_file)
err_psf_image = psf_correction(err_psf_image,setting)

plot_image(err_psf_image,setting,savefig_path+"/err_psf.png",err_image=True)


kwargs_psf = {'psf_type': "PIXEL", 
              'kernel_point_source':psf_image,
              'point_source_supersampling_factor':setting.pssf,
              'psf_error_map':err_psf_image}"""
kwargs_psf = init_kwrg_psf(setting,saveplots=True,backup_path=backup_path)


# ## Modelling

# In[ ]:


# Choice of profiles for the modelling
# Lens (mass) profile with perturber and the Shear
"""lens_model_list = ['SIE']
if CP:
    lens_model_list = ['PEMD']
lens_model_list = [*lens_model_list,'SIS','SHEAR_GAMMA_PSI']

# Lens light profile with perturber
if setting.sub==False:
    light_model_list = ["SERSIC_ELLIPSE", "SERSIC","UNIFORM"]
else:
    light_model_list = ["SERSIC","UNIFORM"]
if hasattr(setting,"no_pert"):
    if setting.no_pert==True:
        light_model_list=["UNIFORM"]
else:
    setting.no_pert=False 
# Source host gal
if not WS:
    source_model_list = ['SERSIC_ELLIPSE']

# Source light model: point 
point_source_list = ['LENSED_POSITION']"""


# In[ ]:


kwargs_params = {'lens_model':        setting.lens_params,
                'point_source_model': setting.ps_params,
                'lens_light_model':   setting.lens_light_params}
if not WS:
    kwargs_params['source_model'] =   setting.source_params


# ### Parameters for the PSO/MCMC runs

# In[ ]:


kwargs_likelihood = {}
kwargs_likelihood["check_matched_source_position"]=True
kwargs_likelihood["source_position_tolerance"]=0.01
kwargs_likelihood["force_no_add_image"] = True 
kwargs_likelihood["check_positive_flux"] = True  


# In[ ]:


# MOD_CUSTOM_LIKE & MOD_LLFR
phi_ll     = setting.phi_ll if setting.sub else None
q_ll       = setting.q_ll   if setting.sub else None
lens_model = init_lens_model(setting)
pycs_magcomb = get_combined_mag(setting.pycs_lensname+"_"+setting.pycs_dataname,main_dir_path=setting.pycs_path)
FR,sig_FR = np.transpose([[get_FR(res),get_sig_FR(res,err)] for res, err in zip(pycs_magcomb.results,pycs_magcomb.tot_error)]).tolist()

Rmag_ABC,sig_Rmag_ABC = [],[]
for i,lbl in enumerate(pycs_magcomb.labels):
    if lbl in ["AB","AC"]:
        Rmag_ABC.append(FR[i])
        sig_Rmag_ABC.append(sig_FR[i])
test_ln  = len(Rmag_ABC) == len(sig_Rmag_ABC)
test_ln2 = len(Rmag_ABC) ==2
if not test_ln or not test_ln2: 
    raise RuntimeError("Something wrong with the Rmag_ABC & sig")
kwargs_likelihood["custom_logL_addition"] = logL_combined(Rmag_ABC,sig_Rmag_ABC,lens_model,SUB=setting.sub,phi_ll=phi_ll,q_ll=q_ll)


# In[ ]:


kwargs_likelihood["image_likelihood_mask_list"] =  [mask.tolist()]


# In[ ]:


"""kwargs_model = {'lens_model_list': lens_model_list,
                'lens_light_model_list': light_model_list,
                'point_source_model_list': point_source_list,
                'additional_images_list': [False], #list of bools (same length as point_source_type_list).
                # If True, search for additional images of the same source is conducted.
                 'fixed_magnification_list': setting.fixed_mag,  # list of bools (same length as point_source_type_list).
                #If True, magnification ratio of point sources is fixed to the one given by the lens model 
                }
if not WS:
    kwargs_model['source_light_model_list'] = source_model_list
"""

######
kwargs_model            = get_kwargs_model(setting)
lens_model_list         = kwargs_model['lens_model_list']
lens_light_model_list   = kwargs_model['lens_light_model_list']
point_source_model_list = kwargs_model['point_source_model_list']
if not WS:
    source_model_list = kwargs_model['source_light_model_list']
######
    
kwargs_numerics = init_kwrg_numerics(setting)
multi_band_list = [[kwargs_data, kwargs_psf, kwargs_numerics]]
# if you have multiple  bands to be modeled simultaneously, you can append them to the mutli_band_list
kwargs_data_joint = {'multi_band_list': multi_band_list, 
                     'multi_band_type': 'multi-linear'  
                     # 'multi-linear': every imaging band has independent solutions of the surface brightness, 
                     #'joint-linear': there is one joint solution of the linear coefficients \
                     # demanded across the bands.
                    }


# In[ ]:


if setting.sub ==False:
    joint_lens_with_light=[[0,0,["center_x","center_y"]],[1,1,["center_x","center_y"]]]
else:
    joint_lens_with_light=[[0,1,["center_x","center_y"]]]


# In[ ]:


kwargs_constraints = {'num_point_source_list': [4], 
                      'solver_type': 'NONE',
                      'joint_lens_with_light':joint_lens_with_light}
# mod free source
if not WS and not setting.FS:
    kwargs_constraints['joint_source_with_point_source'] = [[0, 0]]
                     
#  'joint_lens_with_light': list [[i_light, k_lens, ['param_name1', 'param_name2', ...]], [...], ...],
#   joint parameter between lens model and lens light model


# In[ ]:


if append_MC :
    mcmc_file_name = savemcmc_path+setting_name.replace("settings","mcmc_smpl")+".json"
    with open(mcmc_file_name, 'r') as f:
        mc_init_sample= np.array(json.load(f))
    try:
        mcmc_logL_file_name = savemcmc_path+setting_name.replace("settings","mcmc_logL")+".json"
        with open(mcmc_logL_file_name, 'r') as f:
            mc_init_logL= np.array(json.load(f))
    except:
        mc_init_logL= np.array([])
else:
    mc_init_sample = None
    mc_init_logL= None


# In[ ]:


# Try to solve the "OSError: [Errno 24] Too many open files" by deleting the 
# n_run_cut implementation
from my_lenstronomy.my_fitting_sequence import MyFittingSequence  # ONLY IMPORT IT HERE OR IT BREAK THE CODE

fitting_seq = MyFittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints,\
                              kwargs_likelihood, kwargs_params)
if not append_MC:
    fitting_kwargs_list = [['MY_PSO', {'sigma_scale': 1., 'n_particles': n_particles, 
                                       'n_iterations': n_iterations,"path":savemcmc_path,"threadCount":threadCount}]]
else:
    fitting_kwargs_list = []
    n_burn=0
if RND == False:
    np.random.seed(3)


fitting_kwargs_list.append(['MCMC', {'n_burn': n_burn, 'n_run': n_run, 'walkerRatio': 10, 'sigma_scale': .1,\
                                     "threadCount":threadCount, 'init_samples':mc_init_sample}])
# First chain_list with only the PSO, the burn_in sequence and the first n_run_cut
chain_list = fitting_seq.fit_sequence(fitting_kwargs_list) 
sampler_type, mc_init_sample_i, param_mcmc, mc_init_logL_i  = chain_list[-1]
#append the previous results
if append_MC:
    mc_init_sample = np.array([*mc_init_sample,*mc_init_sample_i])
    mc_init_logL   = np.array([*mc_init_logL,*mc_init_logL_i])
else:
    mc_init_sample = mc_init_sample_i
    mc_init_logL   = mc_init_logL_i

# save here chain_list results
save_mcmc_json(setting=setting,data=mc_init_sample,filename="mcmc_smpl",backup_path=backup_path)
save_mcmc_json(setting=setting,data=mc_init_logL,  filename="mcmc_logL",backup_path=backup_path)
if "PSO" in chain_list[0][0]:
    save_mcmc_json(setting=setting,data=chain_list[0],filename="pso",backup_path=backup_path)
    
kwargs_result   = fitting_seq.best_fit()
param_file_name = savemcmc_path+setting_name.replace("settings","mcmc_prm")+".dat"
with open(param_file_name, 'w+') as param_file:
    for i in range(len(param_mcmc)):
        param_file.writelines(param_mcmc[i]+",\n")


# In[ ]:


# Reconstruct mcmc chain
kw_mcmc        = get_res.get_mcmc(setting_name,backup_path)
chain_list[-1] = ['MCMC',kw_mcmc["mcmc_smpl"],kw_mcmc["mcmc_prm"],kw_mcmc["mcmc_logL"]]


# In[ ]:


print_res.write("kwargs_model:"+str(kwargs_model)+"\n")
print_res.write("kwargs_numerics:"+str(kwargs_numerics)+"\n")
print_res.write("kwargs_constraints:"+str(kwargs_constraints)+"\n")
del kwargs_likelihood["image_likelihood_mask_list"]
print_res.write("kwargs_likelihood:"+str(kwargs_likelihood)+"\n")
if run_type%1==0:
    print_res.write("PSO particles: "+str(n_particles)+"\n")
    print_res.write("PSO run steps: "+str(n_iterations)+"\n")
print_res.write("MCMC run steps: "+str(n_run)+"\n")
print_res.write("MCMC burn in steps: "+str(n_burn)+"\n")
print_res.write("number of non-linear parameters in the MCMC process: "+ str(len(param_mcmc))+"\n")
print_res.write("parameters in order: "+str(param_mcmc)+"\n")
print_res.write("number of evaluations in the MCMC process: "+str(np.shape(chain_list[-1][1])[0])+"\n")
print_res.write("#################################\n")

##Printout of the results with errors
"""
kwargs_sigma_upper={}
kwargs_sigma_lower={}
n_ra=0
n_dec=0
sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]
for i in range(len(param_mcmc)):
    val_min, val, val_max = corner.quantile(samples_mcmc[:,i],q=[0.16, 0.5, 0.84])
    sig_min = np.abs(val_min-val)
    sig_max = val_max - val
    sig_ref = np.min([sig_max,sig_min])
    n_exp = np.log10(1/sig_ref)
    fact=pow(10,int(n_exp)+2)
    if sig_min==sig_max:
        print_res.write(param_mcmc[i]+"  " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " +- "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
    else:
        print_res.write(param_mcmc[i]+" " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " - "+str(np.trunc(np.array(sig_min)*fact)/fact)+\
                        " + "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
    if param_mcmc[i]!="ra_image" and param_mcmc[i]!="dec_image":
        kwargs_sigma_lower[param_mcmc[i]]=sig_min
        kwargs_sigma_upper[param_mcmc[i]]=sig_max

    elif param_mcmc[i]=="ra_image":
        kwargs_sigma_lower["ra_image_"+str(n_ra)]=sig_min
        kwargs_sigma_upper["ra_image_"+str(n_ra)]=sig_max            
        n_ra+=1
    else:
        kwargs_sigma_lower["dec_image_"+str(n_dec)]=sig_min
        kwargs_sigma_upper["dec_image_"+str(n_dec)]=sig_max            
        n_dec+=1
print_res.write("\n#################################\n")
"""
get_res.get_sigma_kw(setting,mcmc_chain=chain_list[-1],print_res=print_res,save=True)


# In[ ]:


# Plot of the obtained models
v_min,v_max     = setting.v_min,setting.v_max
res_min,res_max = setting.res_min,setting.res_max

modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,likelihood_mask_list=[mask.tolist()],\
                      arrow_size=0.02, cmap_string="gist_heat")

if not WS:
    from plotting_tools import plot_model    as PM
else:
    from plotting_tools import plot_model_WS as PM
    
PM(modelPlot,savefig_path,v_min,v_max,res_min,res_max)


# In[ ]:


#Printout of all results after obtaining the amplitude
with open(savefig_path+"/read_results.data","wb") as f:
        pickle.dump(kwargs_result, f)
        
print_res.write("kwargs_results:\n")
for res in kwargs_result:
    len_res = len(kwargs_result[str(res)])
    for i in range(len_res):
            print_res.write(str(res)+" "+str(i)+"\n")
            for j in kwargs_result[str(res)][i]:
                print_res.write(str(j)+": "+str(np.trunc(np.array(kwargs_result[str(res)][i][str(j)])*1e4)/1e4)+"\n")
            print_res.write("\n")
    print_res.write("\n")

print_res.write("\n#################################\n")


# In[ ]:


logL   = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
n_data = modelPlot._imageModel.num_data_evaluate
print_res.write(str(-logL * 2 / n_data)+' reduced X^2 of all evaluated imaging data combined\n')
print_res.write("################################\n")


# In[ ]:


#Normalised plot
f, axes = plt.subplots(figsize=(10,7))
modelPlot.normalized_residual_plot(ax=axes,v_min=res_min, v_max=res_max)
plt.savefig(savefig_path+"normalised_residuals.png")
plt.close()


# In[ ]:


#Caustics
f, axes = plt.subplots(figsize=(10,7))
modelPlot.source_plot(ax=axes, deltaPix_source=0.01, numPix=1000, with_caustics=True)
f.tight_layout()
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
plt.savefig(savefig_path+"caustics.png")
plt.close()


# In[ ]:


f, axes = plt.subplots(figsize=(10,7))
modelPlot.decomposition_plot(ax=axes, text='Point source position', source_add=False, \
                lens_light_add=False, point_source_add=True, v_min=v_min, v_max=v_max)
plt.savefig(savefig_path+"point_source_position.png")
plt.close()


# In[ ]:


#CHECK_FR

#Since the time is the same for all images (considering no time delay, or negligible), we can consider the 
# flux ratio to be amp_i/amp_max

FR,ratio_name = mag_remastered.flux_ratio(setting,kwargs_result,outnames=True)
print_res.write("Flux ratio for "+setting.filter_name+"\n")

for i,FR_i in enumerate(FR):
    print_res.write("Flux ratio:"+str(FR_i)+" "+str(ratio_name[i])+"\n")    
print_res.write("########################\n")


# In[ ]:


# the results of the MCMC chain

i=0
increment= 4
samples_mcmc       = kw_mcmc["mcmc_smpl"]
samples_transposed = np.array(samples_mcmc).T

n_sample = len (samples_mcmc)
samples_mcmc_cut = samples_mcmc
if not len(samples_mcmc) == 0:
    n, num_param = np.shape(samples_mcmc_cut)
    plot = corner.corner(np.array(samples_mcmc_cut[:,:]), labels=param_mcmc[:], show_titles=True)
    plot.savefig(savefig_path+"MCMC_posterior.png")

    #We implement the time delay posterior
    #####################################
    lensModel = init_lens_model(setting)
    #Note: the fermat potential is (should be, double check it) independent on the cosmology, 
    # while the time delay is
    mcmc_new =[]
    labels_new = ["$\phi_{Fermat}A$", "$\phi_{Fermat}C$","$\phi_{Fermat}B$" ,"$\phi_{Fermat}D$","$\Delta t_{A-C}$", "$\Delta t_{A-B}$", "$\Delta t_{A-D}$"]

    for i in range(len(samples_mcmc_cut)):
        #MOD_PROD
        try:
            kwargs_result_i = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
        except NameError:
            print("ERROR: The function setting.produce_kwargs_result is not implemented in ",setting_name,".py file. ")
            raise
        x_image= kwargs_result_i["kwargs_ps"][0]["ra_image"]
        y_image= kwargs_result_i["kwargs_ps"][0]["dec_image"]
        fermat_potential = []
        t_days=[]
        for j in range(len(x_image)):
            fermat_potential.append(lensModel.fermat_potential(x_image[j],y_image[j],kwargs_result_i["kwargs_lens"]))
            t_days.append(lensModel.arrival_time(x_image[j], y_image[j],  kwargs_result_i["kwargs_lens"]))
        mc_i = []
        for j in range(len(fermat_potential)):
            mc_i.append(fermat_potential[j])
        dt=[t_days[0] - t_days[1], t_days[0] - t_days[2], t_days[0] - t_days[3]]
        for j in range(len(dt)):
            mc_i.append(dt[j])
        mcmc_new.append(np.array(mc_i))
    plot = corner.corner(np.array(mcmc_new), labels=labels_new, show_titles=True)
    plot.savefig(savefig_path+"MCMC_time_delay.png")


# In[ ]:


#We save the time delay with the error measure
print_res.write("\n#################################\n")
print_res.write("\nTime delay and fermat potential with MCMC-derived error\n")
kwargs_fermat={}
for i in range(len(labels_new)):
    val_min, val, val_max = corner.quantile(np.array(mcmc_new)[:,i],q=[0.16, 0.5, 0.84])
    sig_min = np.abs(val_min-val)
    sig_max = val_max - val
    sig_ref = np.min([sig_max,sig_min])
    n_exp = np.log10(1/sig_ref)
    fact=pow(10,int(n_exp)+2)
    if sig_min==sig_max:
        print_res.write(labels_new[i]+"  " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " +- "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
    else:
        print_res.write(labels_new[i]+" " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " - "+str(np.trunc(np.array(sig_min)*fact)/fact)+\
                        " + "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
    kwargs_fermat[labels_new[i]]=(val,sig_min,sig_max)
print_res.write("\n#################################\n")


# In[ ]:


#Time delay
print_res.write("\n########################\n")
print_res.write("Time delay\nWe are considering the standard cosmology of astropy\n")

# if we consider the new image we convert the pix_coord of the objects
# and after that obtain the angle via pix_scale

lensModel = init_lens_model(setting)

# here are the predicted time delays in units of days
x_image= kwargs_result["kwargs_ps"][0]["ra_image"]
y_image= kwargs_result["kwargs_ps"][0]["dec_image"]
t_days = lensModel.arrival_time(x_image, y_image, kwargs_result["kwargs_lens"])
#This are "absolute time", useless since we cannot measure them

# here are the predicted time delays in units of days
t_days =[]
k=0
for x in x_image:
    t_days.append(lensModel.arrival_time(x, y_image[k],  kwargs_result["kwargs_lens"]))
    k+=1
    
print_res.write ('time delays: '+ str(t_days[0] - t_days[1])+" "+str(t_days[0] - t_days[2])+"  "+\
                  str(t_days[0] - t_days[3])+"\n")

print_res.write ("Between images 0-1, 0-2, 0-3 with coord x:"+str(x_image)+" \nand coord y"+str(y_image)+"\n")


# In[ ]:


#MOD_SOURCE 
from source_pos import get_source_pos_MCMC
kwargs_source,str_src = get_source_pos_MCMC(setting,svfg=True)
print_res.write(str_src)


# In[ ]:


#Closing result.txt
print_res.close()


# In[ ]:


# I save the kwargs result in a pickly, readable way
def pickle_results(res,name):
    if not ".data" in name:
        name+=".data"
    with open(savefig_path+name,"wb") as f:
        pickle.dump(res, f)
last_kw = {"read_results":kwargs_result,
#           "read_sigma_up":kwargs_sigma_upper,
#           "read_sigma_low":kwargs_sigma_lower,
           "read_fermat":kwargs_fermat,
           "read_source":kwargs_source,
           "FR":FR}
for nm in last_kw:
    pickle_results(last_kw[nm],nm)
    


# In[1]:


# add some final commands

import last_commands
    
for i in last_commands.progs:
    last_commands.last_command(setting_name, i,log=True) 
    
    
from check_success import check_success
check_success(setting_name,verbose=1)


# In[ ]:


success(sys.argv[0])

