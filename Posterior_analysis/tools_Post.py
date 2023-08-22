import copy
import numpy as np

def get_array_BC(array,relation="ratio"):
    # assuming array ordered in AB,AC,AD
    if relation=="ratio":
        array_cp    = np.array(copy.deepcopy(array))
        # AC/AB = (C/A)/(B/A) = (C*A)/(B*A) = C/B = BC
        array_BC    = array_cp[1]/array_cp[0]  
        array_cp[2] = array_BC 
        array       = array_cp #AB,AC,BC
    elif relation=="subtraction": 
        array_cp    = np.array(copy.deepcopy(array))
        # AC-AB = (C-A)-(B-A) = C - A - B + A = C - B = BC
        array_BC    = array_cp[1] - array_cp[0]  
        array_cp[2] = array_BC 
        array       = array_cp #AB,AC,BC
    else:
        raise ValueError(f"Argument relation should be either 'ratio' or 'subtraction', not {relation} ")
    return array