import numpy as np
from MtL_ratio import Grid_Class
import matplotlib.pyplot as plt

#precisions = [0.05,0.1,0.5,0.7,0.9]

precisions = [0.05,0.06,0.07,0.08,0.09,0.1]#,0.5,0.7,0.9]
radius = np.arange(1,20,.5)
"""

for prec in precisions:
    A_p = []
    for r in radius:
        gr = Grid_Class(10,10,r,prec)
        _,area_circ = gr.circularise_flat_grid_with_area()   
        A_rp = np.sum(area_circ)
        a_real = np.pi*(r**2)
        print("r =",np.round(r,2),"A_th = ",np.round(a_real,2),"A_meas = ",np.round(A_rp,2)," ratio meas/real = ",np.round(A_rp/a_real,2)," prec =",prec)
        A_p.append(A_rp)
    plt.scatter(radius,A_p,label=str(prec),marker=".")
plt.legend()
plt.savefig("test/precision_MtLratio.png")
plt.close()

"""

import timeit
for prec in precisions:
    A_p = []
    start = timeit.default_timer()
    for r in radius:
        gr = Grid_Class(10,10,r,prec)
        _,area_circ = gr.circularise_flat_grid_with_area()   
        A_rp = np.sum(area_circ)
        a_real = np.pi*(r**2)
        print("r =",np.round(r,2),"A_th = ",np.round(a_real,2),"A_meas = ",np.round(A_rp,2)," ratio meas/real = ",np.round(A_rp/a_real,2)," prec =",prec)
        A_p.append((A_rp-a_real)/a_real)
    stop = timeit.default_timer()
    plt.scatter(radius,A_p,label="<prec "+str(prec)+">: "+str(abs(np.round(np.mean(A_p),4))*100)+"% in "+str(np.round(stop-start,3))+" sec",marker=".")
    
plt.legend()
plt.ylabel(" A_meas - A_true / A_true") 
plt.savefig("test/precision_MtLratio_.png")
plt.close()

