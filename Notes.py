#!/usr/bin/env python
# coding: utf-8

# ## 2th January 2023
# - so rn 1 big problem and 1 minor
# - big problem:
#     - the f105w CPRP is pretty much off
#     - ![image-4.png](attachment:image-4.png)
#     - ![image.png](attachment:image.png)
#     - having checked the difference between that and all the rest, we found even more disagreements
#     - also, the results obtained now with the same setting file (and a longer run) are slightly different:
#         - ![image-5.png](attachment:image-5.png)
#         - possibly the problem is connected to not obviously not converged mcmc runs
#         - or the Main program had slightly changed since last time I run it and now it's a different result:
#         - ![image-3.png](attachment:image-3.png)
#         - still they both are far from the CPRP
#         - ![image-6.png](attachment:image-6.png)
#         
#     - the point can be that we don't need/have to describe it as CPRP does, ie not EXACTLY the same prior, as the requirements for the precision and the light profiles can be different
#         - also, are we sure that this prior doesn't get affected by the orientation of the image? 
#         - e1 and e2 mainly for the image
#         - no not anymore, as they are corrected for the CDij matrix right
#     - so redo ALL with same LENGHT of MCMC and see what happens
#     - checked the "genealogy" of setting to see where and which choices were done (see sett_genealogy_tree.py)
#     
# - minor problem:
#     - we have a "hole" in the observation of 2 month:
#     - ![image-2.png](attachment:image-2.png)
#     - 1 month we couldn't observe, the second all observations were "bad"
#     - check this, it's suspicious
# 
# 
# - create program to print a latex table with the prior values 

# ## 3th January
# - while writing about the Prior i realised something:
#     - the source pos is given in the prior, but actually it doesn't matter bc it is joint with the point source
#  >
#     hierarchy is as follows:
#     1. Point source parameters are inferred
#     2. Lens light joint parameters are set
#     3. Lens model joint constraints are set
#     4. Lens model solver is applied
#     5. Joint source and point source is applied
# - so it looks like (check also Lenstronomy/Sampling/parameters.py ) it looks like it's ignored and only inferred from the backward ray tracing equation for the images position (unless FS=True)
# - was redoing the plot for "Light Profiles and Lens Light Model", with the lens ligh model for f814w:
#     - ![image.png](attachment:image.png)
#     - the red line is weighted average, the green is average
#     - the two vertical lines approx indicate the distance from the center for the closest and furthest image (ie the region we study in the lens modelling)
#     - about the ellipticity, what we have is:
# >     def logL_ellipt_q(self,q,q_ll,sigma_q=.1):
# >        ## follow Note 27th June
# >        diff = q-(q_ll-0.1)
# >        if diff>=0:
# >           return 0
# >        else:
# >            return -diff**2/sigma_q**2/2
# - this is due to the constrain  $𝑞 \geq 𝑞_L − 0.1$. and given that the sigma is .1, it still cover most of the uncertainty of eps. so we shouldn't worry too much
#     - rather, are we sure is in anyway useful, that it does really change something?
#     - from my experience, it doesn't
#     
# - tomorrow: re-read and finish at least 2 parts and send them to stella &Arno
# 
# ## 4th Jannuary
# - read "Mass-sheet degeneracy: Fundamental limit on the cluster mass reconstruction from statistical (weak) lensing" partially
# - to read "The probability of galaxy-galaxy strong lensing events in hydrodynamical simulations of galaxy clusters"
# - corrected 1st section - ready to send
#     - add 1 color image
#     
# ## 5th Jannuary
# - asked Raphael about how to produce nice color images
#     - actually would be cool to have AI do that
# - check results
# - read meneghetti
# - write 2 sections
# - since the f140w is actually slightly changed in the new result, In order to avoid problems I saved the previous results as "mcmc\_...\_old{\_ws}" to check how much/if it changes something
#     - for f814w_CP_logL2_IV_ws  f160w_7030_CP_logL2_mkSource_IV_SC  f475w_CP_logL2_ws  f160w_CP_logL2_mk32 
# - for now the result from f140 are better in agreement with f105 (the longer run version) then before
# 
# ## 9th Jannuary
# - ![image-2.png](attachment:image-2.png)
# - they are in slight to severe tension
# - mainly f160w_7030_CP_logL2_mkSource_IV_SC
# - but overall, it looks like the infrared got pushed to the bottom left, while the optical did not change much
# - why:
#     - something changed on the code?
#         - maybe ellipticity?
#             - should not have changed in a while
#     - longer run reached better convergence?
#         - this would be bad
#         - worse even:
#             - log_f140w_CP_logL2_SC_I_2Jan.dat:Max iteration reached! Stopping.
#             - log_f160w_CP_logL2_mk32_2Jan.dat:Max iteration reached! Stopping.
#             - log_f105w_CP_logL2_skycorr_tmp_ws_23Dec.dat:Max iteration reached! Stopping.
#         - while:
#             - log_f160w_7030_CP_logL2_mkSource_IV_SC_2Jan.dat:Converged after 796 iterations!
#             - log_f475w_CP_logL2_ws_2Jan.dat:Converged after 364 iterations!
#             - log_f814w_CP_logL2_IV_ws_2Jan.dat:Converged after 377 iterations!
#             - 7030 seems very similar to its version CPRP
#             - ![image-3.png](attachment:image-3.png)
#             
# ## 10th Jannuary
# - mcmc_lenght = mcmc_run_steps\*run_factor\*n_params\*walker_ratio ?
# - but not clear why there is still a factor of 10 ...:
# - f160w_7030_CP_logL2_mkSource_IV_SC:
#     - lenght = 4320000
#     - mcmc_run_steps = 800
#     - rf = 2 (default)
#     - n_params = 27
#     - walker_ratio = 10
# - somewhere, a factor of 10 appears...
# - difference of f160w's:
#     - already done previously
#     - repeated with a small change (before I was changing the header description and cutting, now I used radeccuter):
#     - the difference are present, especially at the image pos (can be ignored)
#         - should be deconvolved by the psf to be actually comparable, but that would be complicated
#         
#     - and more interestingly at the source light pos (weakly more light in the 
#     - <img src="notes_images/diff_f160w.png"  width="550" height="400">
#     - this is $ (F160W - fact *F160W_{7030})/ \sqrt{Err_{F160W}^2+ (Err^{7030}_{F160W}*fact)^2} $
#     - so F160w
# - found that they implemented mcmc Zeus in lenstronomy -> faster? 
#     - to be tested, but not now
#         - check /home/jackquei/Desktop/UNI_Muenchen/Thesis/lenstronomy/Notebooks/test_zeus.ipynb
#         - looks like it's slightly better but slower
#         - not good
#     - only if we obtain good results (and identical to emcee) faster 
#     

# ## 11th Jan
# - write
#     - correct elena's point
#     - write/correct 1 chap: lens modelling
# ## 12th Jan
# - limited by battery:
#     - stay @office for the call and then leave
#     - or go back for the call
# - f160w_7030:
#     - ![image.png](attachment:image.png)
#     - ![image-3.png](attachment:image-3.png)
#     - the only damn difference is that the blue one has been running 2 times longer
# - f814w: 
#     - reach convergence
#     - but q is very much at the boundary of the allowed parametric space
#     - ![image-2.png](attachment:image-2.png)
# - rubber ducking
# 
# - a lot of points:
# <img src="notes_images/120123.jpeg"  width="550" height="400">
#     - main point: some model works (are constant with number of mcmc steps) some other aren't (f814 and f160w_7030) which is odd : the model might be correct, but the data is messed up?
#     - checked the psf for f160w_7030:
#     
# ![image-4.png](attachment:image-4.png)
# - psf of 7030 is not centered!!
# - used the following (ugly) code:

# In[ ]:


data = kw_psf_7030['kernel_point_source']  ;  cnt_data_x = int(len(data[1])/2.) ; cnt_data_y =  int(len(data[0])/2.) ; dx  = int(cnt_data_x/2) ; x0  = int(cnt_data_x-dx) ; x1  = int(cnt_data_x+dx) ; proj_data  = data[cnt_data_y][x0:x1] ; x_projdata = np.arange(x0,x1,dtype=int) ; plt.plot(x_projdata,proj_data,c="r",label="f160w_7030, cnt:"+str(cnt_data_x)) ; plt.axvline(cnt_data_x,c="r",linestyle="--") ; data = kw_psf_160['kernel_point_source']  ;  cnt_data_x = int(len(data[1])/2.) ; cnt_data_y =  int(len(data[0])/2.) ; dx  = int(cnt_data_x/2) ; x0  = int(cnt_data_x-dx) ; x1  = int(cnt_data_x+dx) ;  proj_data  = data[cnt_data_y][x0:x1] ; x_projdata = np.arange(x0,x1,dtype=int) ; plt.plot(x_projdata,proj_data,c="b",label="f160w, cnt:"+str(cnt_data_x)) ; plt.axvline(cnt_data_x,c="b",linestyle="--") ; plt.legend() ; plt.savefig("test/psf_f160w_superposed.pdf") ; plt.close()


# - check f814w vs f475
#     - standard:
#     - ![image.png](attachment:image.png)
#     - normalised by the total summ :
#     - ![image-3.png](attachment:image-3.png)
# - a small difference of the wings, but for the rest they are very similar
# - maybe the lens light subtraction went wrong?
# 
# 
# - try f160w psf for 7030
# - investigate PSO
#     - pso with different params for convergence
#     - check if different conv params affect the result/the steps
#     - 20 different runs -> check if params ends up being the same and are convergent
#         - with different random seeds, all the rest is the same
#         
# - rewrite mcmc_behaviour_plot.py better
# 

# ## 17th Jannuary
# - SGL of CMB?
# - check lens light models for f814 and f475 and compare with lens light models from f160 and others
#     - ![image.png](attachment:image.png)
#     - rel. well in agreement (left 814, right f475)
#     - check with f160w, f160w_7030, f105 and f140
#######################################
f160w_7030_CP_logL2_mkSource_III_SC
#######################################
e1 0.2313287732613967
e2 -0.07119780673041122
mcmc res
e1 0.13827307957154733
e2 -0.018326985382265782
#######################
f160w_CP_logL2_mk32
#######################
e1 0.23836997162110912
e2 -0.06562278360170291
mcmc res
e1 0.160726053255367
e2 -0.023834127203254347
#######################
f140w_CP_logL2_SC_I
#######################
e1 0.23757034387657303
e2 -0.0674878376700991
mcmc res
e1 0.17442800190271526
e2 -0.046585708534047206
#############################
f105w_CP_logL2_skycorr_ws
#############################
e1 0.19934089340800729
e2 -0.05620205473154379
mcmc res
e1 0.15993500458844767
e2 -0.06473144037636824

# ## 18th Jan 
# - PSO reach convergences in different points:
#     - pso_it_0: conv after 313 it
#     - pso_it_1 : conv after ~330 it
# - check what is the value obtained from the kwargs_result and where it appears in the pso
#     - then print out with pso_check_res these values for the various pso chains 
#     - and try to find out why they are different
#     
# - ABOUT THE PSO:
#     - if there is only that, the kwargs_result is == best_fit == global_best
#     - this should be the last step (or not too far away from it?)
#         - testing
#         - so:
#             - for 2 out of 3 (and the 3rd was an outlier), the result outputed was = to global best (=min -log_likelihood)= to the latest step:
#             - have to plot the "directions" of the steps
#             
#     - idea to visualise how it works for the f814:
#         - select one "reference point" in the parametric space
#                 - result from previous run?
#                 - other results from the other filters?
#         - plot the modulus of the "distance" from it over the different steps
#         - ![image-2.png](attachment:image-2.png)
#         - for now only 2 pso runs, but both are interestingly going far away from the reference point
#         - both are converged
# - following cell: analysis of the lightcurve braking
- find out why lightcurves brakes :
    - from 20220428 to 20220618
        - 3 obs day
           - 2 times our was observed 0429, 0501
           - both were reduced
                - 220501_091 220501_092 220501_093 are good enough
                    - one final stacked image present
                    - transp: 0.703 (min = .7 -> fine)
                    - seeing: 0.976 (
                    - between a set of 12 images that do not end up in the final list
                    - 20201115 20210503 20211010 20211110 20211209 20220223 20220224 20220228 20220501 20221104 20221114 20221121
                    - because of this error:
dividing flxscale = 
Traceback (most recent call last):
  File "/data/wst/u/wst/WWFI/SCRIPTS/subtractfits.py", line 83, in <module>
    mdata /= mhdr['FLXSCALE']
TypeError: ufunc 'true_divide' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''
mv: cannot stat ‘seSDSSJ1433_V200528_g_156.fits’: No such file or directory
mv: cannot stat ‘e.seSDSSJ1433_V200528_g_156.fits’: No such file or directory
                    - the fluxscale should be given by the qdiffima3.sh program
                        - it reads it from the SDSSJ1433_V200528_g_$NUM*_newzp.txt file
                            - absent! -> absent for all the 12 images
                                - produced by the /data/wst/u/mkluge/scripts/psfphotometry_wrapper.py program
                                - log files all have this error "ValueError: Illegal value: masked." -> all images have this error
                                - but they also have at least some that ended well
                                - actually the psf is produced for all of them
                                - and the psfphotometry logs seems quite similar to the standard ones.... apart the fact that in average there are more log files then the normal
                                -                                 
    - from 20220821 to 20221102
        - 39 obs days
            - 31 with SDSS observations
            - 10 still to reduce = 20220823 20220824 20220825 20220828 20220829 20220903 20220904 20220908 20220922 20220923
                - running
            - 21 had their directory -> none was fully reduced
                -  20220901 20220910 20220911 20220912 20221003 20221005 20221006 20221007 20221009 20221015 20221016 20221017 20221019 20221022 20221026 20221027 20221028 20221029 20221030 20221031 20221101
                - none were stacked
                - check if mask was done 
                - actually they had problems:
                    - log.autoreduce_V2_221102_194812: 
> 251 jobs done in 39 minutes.

ERROR: 1 out of 251 frames have not been reduced.

# ## 20th Jannuary
# - check dates
# - pso
#     - continue running for f814w
#         - actually it is still running...
#     - replotting with 3 psos
#     - consider to plot the pos of only the parameter e1 (or was it e2) over time wrt 
#         - pso
#         - mcmc
# - to complete: 
#     - implement averaged_plot in all behaviour plots
# ## 23th Jannuary
# - usm machine not reacheable
#     - now yes
#     - the f160w_7030_CP_logL2_mkSC_psf_long is actually more similar to f160w_7030_CP_logL2_mkSource_IV_SC then f160w_7030_CP_logL2_mkSC_psf
#         - ![image.png](attachment:image.png)
#         - as expected the slight modification of the PSF did not change the result
#     - behaviour of e1 for f160w_7030_psf_long
#     - ![image-2.png](attachment:image-2.png)
#     - left old, right new:
#     - ![image-3.png](attachment:image-3.png)
#     - and the old e1 behaviour:
#     - ![image-5.png](attachment:image-5.png)
#     - not converged
#     - comparison btw bad (top f160w_psf_long) and good(bottom f475w_CP_long) e1 bhv
#     - ![image-4.png](attachment:image-4.png)
#     - f160 have a semistable solution at steps~400 which is somewhat close to the actual result of f475. then it leaves it and reach a much lower solution, around which the mcmc has to stay. 
#         - the pso varies on ~1 order of mag larger then mcmc (0.15 vs 0.015)
#             - FOR THIS REASON the convergence of PSO is paramount, else the MCMC will sampled always new spots
#             
#         - depends also on the given parameters for the runs
#         - note: f160w is not converged, while f475 is 
#     - produced images at a given step of the PSO, but hard to see, double check
# - not reduced nights: 
#     - might have been due to lacking of masterbias for that month (now present)
#     - rerunning
# 

# ## 26th Jan
# - for some reason, it starts the PSO with the wrong value for e1:
# - ![image.png](attachment:image.png)
# - not only f814w but all of them
# - and not only for e1 (although e1 is the worts one, as it goes all the way to the upper bound)
#     - in theory (I tested with a very short test) the first step of the PSO is the initial value
#     - but when taking the first step of the actual analysis, they are off from the initial value, for e1 by a lot (reaching maximum allowed value)
#     - it looks like they actually do some steps before, while in the test they don't
#     - why/how? and why does it reach the max value (for f814w and f160w_7030, but also other vary a lot)
# - so, first thing:
#    - lower and upper limit are not inherited directly, but estimated from the sigma given by the kwargs
#        - changed in the newest version of my_fitting_sequence 
#            - could it be a problem doing it this way?
#            - the one given are, by construction, at the same distance (sigma) from the initial position, meaning that the recovered initial position is identical to the given one -> cannot be that
#            - MOD_PSOLIM: corrected them to be the one given by me
#    - the initial value is, through the updateManager, inherited directly. So there is no way (in theory) that it get's out wrong
#        - how can this happen and where
#            - even if initialising the PSO with many more particles, the initial pos is still the one given at the beginning
# - question: PSO by construction have many particles, how are those connected to the single path that we get in the results?
#     - simply they are in a for loop-> this might be why we get a non identical initial pos:
#         - because it is the position of a random particle, which is obtained by sampling
#         - from my_pso.py:
    def sample
        pos=[]
        lkl=[]
        while True:
            for particle in self.swarm:
                pos.append(np.array(particle.position))
# - at the same time:
#     - why exactly at the boundary
#     - why the mean seems to be off -> check that
#     - why when doing the few steps test (even when stopping)
#         - actually this is weird, the conversion gave me exacty the initial value,
#         - and the saved version once gave me another weird point (e1~.4) and another time the initial value
#         - repeated 3 times, 1 time it was different, 2 times it was the same as the initial value
#         - so it might be that it is random...? Sometime it's the same, sometime it isn't
#         - but why when it wasn't, was it at the boundary?
# - actually another question:
#     - since the steps should be the steps of every particle, how come that the output is only as long as the number of steps, and not stepsxn_particles?
#         - ACTUALLY the output of optimise is the following:

# In[1]:


chi2_list = []
vel_list = []
pos_list = []

num_iter = 0
for _ in self.sample(max_iter, c1, c2, p, m, n, early_stop_tolerance):
    chi2_list.append(self.global_best.fitness * 2)
    vel_list.append(self.global_best.velocity)
    pos_list.append(self.global_best.position)
    num_iter += 1
    
    if verbose and self.is_master():
        if num_iter % 10 == 0:
            print(num_iter)

return self.global_best.position, [chi2_list, pos_list, vel_list]


# - meaning that for every step of the sample, only the global best is considered 
#     - initially this is the initial point, unless by chance a single point ends up having a higher logL then the initial
#     - but then how can it be that the first local optima is right at the boundary? 
#         - something is off here
#         - possible: the limits are independent of the actual limits, and the pso sample outside of the limits, for some reason giving the best pos somewhere there?
#             - a) sketchy
#             - b) even with the corrected version I obtained something like that
# - note: while the PSO by itself is not constrained to stay within the boundaries, it call the Likelihood function, which in turn give a huge penalty if it is outside the boundaries
#     - it is not the case here, BUT this could be a problem if the PSO computed boundaries are bad and the particle spawn all outside of them: then they all have the same likelihood.
#         - not this case, also it should be quite easy to see?
#         - as soon as one particle enters back in the region, it should draw all of the others with her
#         - asked in the slack channel
#         - ![image.png](attachment:image.png)
#         - so we have implemented the program in a improper way?
#         - three parameter have the problem of the "sigma boundary" (initial value +- sigma) being "larger" than the forced boundary: theta_E, e1,e2
#             - by quite a lot:
#             ![image-2.png](attachment:image-2.png)
# - do a PSO dummy test:
#     - 1 single param, in a simulated case 
#     - see what happens with the sigma being too large
# - actually in that case probably it woudln't be a problem:
#     - in theory my approach should be sounder against observer bias, but in reality it does spread too much the PSO particles -> meaning that if one finds a good enough local minimum, it will get stuck there and attract the other particles, which might not find a larger minima before getting there
#     - dummy test with many params 
#     - run, it looks like in average it doesn't change much
#     - test f814w with this approach
#         - corrected theta_E, e1, e2
#         - same_prior_bound.py

# ### 1st February
# - make the stock of the situation
#     - to do:
#         - git everything
#         - checked the difference between the light produced by the model for the opticals and the infrared:
#             - actually doesn't tell much... by eye the might look somewhat compatible, also considerign that the Rsersic and the n_sersic are placeholders for the optical as we only have isophotes
# - big big mistake found:
#     - the ellipticity in matthias "isophotes.py" might (probably) be saved as 1-b/a, where a and b are major and minor axis, respectively
#     - but i was reading it as b/a directly
#     - the difference was that the logL term to controll q was completely useless:
#         - logL = 0 if Dq=q-(qll-0.1)>0 -> q>0.5, q_ll (wrong) = 0.3 -> Dq>0
#         - now qll ~ 0.7, so q must be higher then ~0.7 -> which would be more in agreement with the other infrered filters (both f814w and f475w)
#         - still, might be less strong then other problems
#         - note! corrected in the code (extract_q_ll), but in the meantime f814w_CP_logL2_VII_ws has run with the correction defined in the setting file, be aware!
# 
# 
# ### 2nd February
# - little done
# - to do:
#     - create better structure for my code
#     - create definitive prior and relaunch better done f475w and f814w
#         - check redone prior: either that or a similar structure
#     - read/study something
#         - ml? article?
