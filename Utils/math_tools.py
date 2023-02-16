import numpy as np
def sqrt_sum(a,b):
    if np.shape(a)==():
        return np.sqrt(a**2 + b**2)
    else:
        if len(a)!=len(b):
            raise ValueError("Len of a and b must be equal, not: "+str(len(a))+" and "+str(len(b)))
        return np.array([np.sqrt(a[i]**2 + b[i]**2) for i in range(len(a))])
    
def sqrt_sum_list(list_ab):
    if all([np.shape(ab)==() for ab in list_ab]):
        return np.sqrt(np.sum([ab**2 for ab in list_ab]))
    else:
        if any([len(list_ab[0])!=len(ab) for ab in list_ab]):
            raise ValueError("Len of each element must be equal")
        return np.array([np.sqrt(np.sum([ab[i]**2 for ab in list_ab])) for i in range(len(list_ab[0]))])
