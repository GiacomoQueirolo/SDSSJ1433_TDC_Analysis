#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# A general way to order the images in a standardised order


# In[38]:


import numpy as np
#wo A order   x    B       ,     C      ,    D
ra_ordered  = [-0.00221745, -0.76000637, 2.04098104 ] 
dec_ordered = [-3.75414086, -2.12910109, -2.17638207 ]

def image_order(ra,dec,ret_order=True,verbose=True):
    if len(ra)!=len(dec):
        raise ValueError("Give array of the same lenght!")
    if len(ra)==4:
        if verbose:
            print("Consider the first image to be A by construction")
        ra = ra[1:]
        dec= dec[1:]
    image_order = []
    for r,d in zip(ra,dec):
        dist =[]
        for ro,do in zip(ra_ordered,dec_ordered):
            dist.append(np.sqrt((r-ro)**2 + (d-do)**2))
        #print(r,d,dist,np.argsort(dist))
        #min_dist = np.argsort(dist)[0]
        min_dist = np.where(dist==min(dist))[0][0]
        image_order.append(min_dist)
    if not np.all(np.unique(image_order,return_counts=True)[1]==np.ones_like(image_order)):
        raise ValueError("Something went wrong, there are multiply images for the same name")
    
    image_order= np.argsort(image_order)
    
    if ret_order: # I only return the order to be applied to the array 
        return image_order
    else:
        return ra[image_order],dec[image_order]


# In[ ]:


from Utils.tools import get_setting_module
from Utils.get_res import get_mcmc_prm,get_mcmc_smpl

def get_new_image_order(setting,mcmc=None,starting_from_A=True,backup_path="backup_results"):
    setting = get_setting_module(setting,1)
    param_mcmc   = get_mcmc_prm(setting,backup_path)
    if mcmc is None:
        mcmc = get_mcmc_smpl(setting,backup_path)
    len_radec = min([300,len(mcmc)])
    if len_radec<20:
        print("Warning: You are only using "+str(len_radec)+" steps to determine the order.\
                     I hope you know what you are doing")
    ra_ord,dec_ord = [],[]    
    for i in range(len_radec):
        kwargs_result_i = setting.produce_kwargs_result(mcmc,param_mcmc,i)
        x_image         = kwargs_result_i["kwargs_ps"][0]["ra_image"]
        y_image         = kwargs_result_i["kwargs_ps"][0]["dec_image"]
        ra_ord.append(x_image)
        dec_ord.append(y_image)
    ra_ord  = np.mean(ra_ord,axis=0)
    dec_ord = np.mean(dec_ord,axis=0)
    new_order = image_order(ra_ord,dec_ord)
    if starting_from_A:
        new_order+= 1 #bc it gives the order assuming A in 0
        new_order = [0,*new_order]
    return new_order


# In[45]:


# test:
#      A,   C , D  ,  B
"""
ra_names = np.array(["A","C","D","B"])
ra  = np.array([0,-0.79,2.0,  0.01  ])
dec = np.array([0,-2.1,-2.4,-3.8 ])
order = image_order(ra,dec)

order+=1
order= [0,*order]
ra_names[order]
"""

