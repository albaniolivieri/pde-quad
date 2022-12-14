import signal
import time
import math
from sympy import *
from .var_selection import prop_new_vars
from .check_quad import get_quad

class PolynomialSystem: 
    
    def __init__(self, dics, poly_syms, order, var_indep, pde_eq):
        self.dic_t = dics[0]
        self.dic_x = dics[1]
        self.poly_vars = poly_syms
        self.order = order
        self.var_indep = var_indep
        self.pde_eq = pde_eq

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

# idea: structure that has all the information of poly system 
def bnb(new_vars, best_nvars, poly_syst):
    if len(new_vars) >= best_nvars:
        return None, math.inf, 1
    new_vars_named = [(symbols(f'w_{i}'), pol) for i, pol in enumerate(new_vars)]
    result_quad = get_quad(poly_syst.dic_t, poly_syst.dic_x, new_vars_named, poly_syst.pde_eq, 
                           poly_syst.order, poly_syst.var_indep, poly_syst.poly_vars)
    if result_quad[0]:
        return new_vars, len(new_vars), 1
    
    min_nvars = best_nvars
    best_quad_vars = None
    
    traversed_total = 1
    prop_vars = prop_new_vars(result_quad[1], new_vars)
    #print('prop_vars', prop_vars)
    for p_vars in prop_vars: 
        #print('p_vars', p_vars)
        print('new vars added', new_vars + list(p_vars))
        quad_vars, nvars, traversed = bnb(new_vars + list(p_vars), min_nvars, poly_syst)
        traversed_total += traversed
        print('nodes traversed:', traversed_total)
        if nvars < min_nvars:
            min_nvars = nvars
            print('Best quadratization until now:', min_nvars, quad_vars)
            best_quad_vars = quad_vars
    
    return best_quad_vars, min_nvars, traversed_total
            
            