import signal
import time
import math
from sympy import *
from .PolySys import *
from .var_selection import *
from .utils import powerset

ALGORITHM_INTERRUPTED = False

# def signal_handler(sig_num, frame):
#     global ALGORITHM_INTERRUPTED
#     print("The algorithm has been interrupted. Returning the current best.")
#     ALGORITHM_INTERRUPTED = True

# signal.signal(signal.SIGINT, signal_handler)    

def pruning_rule_nvars(nvars, global_nvars):
    if nvars >= global_nvars: return True 
    return False

def pruning_rule_time(start_time, max_time):
    if time.time() - start_time > max_time: return True
    return False

def shrink_quad(quad_vars, poly_syst):
    final_vars = quad_vars
    subsets = powerset(quad_vars)
    for var_group in subsets: 
        poly_syst.set_new_vars(var_group)
        res, _ = poly_syst.get_quad() 
        if res:
            return list(var_group)
    return final_vars

def bnb(new_vars, best_nvars, poly_syst, sort_fun):
    if len(new_vars) >= best_nvars:
        return None, math.inf, 1
    
    poly_syst.set_new_vars(new_vars)
    result_quad = poly_syst.get_quad()
    
    if result_quad[0]:
        shrinked_quad = shrink_quad(new_vars, poly_syst)
        return shrinked_quad, len(shrinked_quad), 1
    
    min_nvars = best_nvars
    best_quad_vars = None
    traversed_total = 1
    prop_vars = prop_new_vars(result_quad[1], new_vars, sort_fun)
    
    for p_vars in prop_vars: 
        quad_vars, nvars, traversed = bnb(new_vars + list(p_vars), min_nvars, poly_syst, sort_fun)
        traversed_total += traversed
        if nvars < min_nvars:
            min_nvars = nvars
            print('Best quadratization until now:', min_nvars, quad_vars)
            best_quad_vars = quad_vars
    
    return best_quad_vars, min_nvars, traversed_total