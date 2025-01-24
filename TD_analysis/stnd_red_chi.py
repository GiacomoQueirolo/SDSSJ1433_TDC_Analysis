import numpy as np

def red_chi(data,sigma,model,n_int_knt,mlparams,n_lc):
    #CHI
    ####
    chi = np.sum(((data - model)/sigma)**2 )
    ####
    
    #DOF
    ####
    #- n_int_knt is the number of INTERNAL knots, the only that can vary
    #- intr. spline always deg = 3 -> meaning f(x) = a +b*x + c*x² + d*x³ -> 4 DoF
    #- for every spline, the n* of DoF is given by: 4 param at the first interval
    # + 1 free param for every new interval (due to constrains at the internal knots)
    # - the internal knot position (both intr and ml) are also free to vary -> int_knot are also DoF
    # - n_intervals = n_int_knt + 1
    
    #Intrinsic spline
    ##
    free_par_spl = 4
    n_free_par = free_par_spl +n_int_knt  #spline params= 4(1st interval) + 1*new_intervals(=n_int_knots)
    n_free_par += n_int_knt # internal knots
    ##
    
    # Mag and time shift for every lc -1
    n_free_par += 2*(n_lc-1)
    ###   
    
    #ML
    ##
    mllist_name = mlparams["mllist_name"]
    if mllist_name is not None:
        for i in mllist_name:
            if mlparams["mltype"]=="polyml":
                n_free_par += mlparams["mlfp"] # free params of the ml polynomial 
            elif mlparams["mltype"]=="splml":
                n_free_par += free_par_spl 
                n_free_par += mlparams["nmlspl"] #n* int knots
                if mlparams["forcen"] is False:
                    n_free_par += mlparams["nmlspl"]
    ##

    DoF = len(data) - n_free_par
    
    #CHI_RED
    ###
    chi_red = chi/DoF
    ###
    return chi_red 

def get_chi_red(spline_i,mlparams,nl_lc):
    spl_mask    = spline_i.datapoints.mask
    m, merr     = spline_i.datapoints.mags[spl_mask], spline_i.datapoints.magerrs[spl_mask]
    f           = spline_i.eval(spline_i.datapoints.jds[spl_mask])
    n_int_knt   = spline_i.getnint()-1 # n* intervals (of intrinsic spline)
    spl_chi_red = red_chi(m,merr,f,n_int_knt,mlparams, nl_lc)
    return spl_chi_red

