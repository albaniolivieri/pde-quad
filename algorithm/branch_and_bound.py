import signal
import time
import math
from .var_selection import prop_new_vars
from .utils import get_diff_order, shrink_quad

ALGORITHM_INTERRUPTED = False

# def signal_handler(sig_num, frame):
#     global ALGORITHM_INTERRUPTED
#     print("The algorithm has been interrupted. Returning the current best.")
#     ALGORITHM_INTERRUPTED = True

# signal.signal(signal.SIGINT, signal_handler)    

def pruning_rule_nvars(nvars, global_nvars):
    """Pruning rule based on the number of variables in the quadratization found.
    
    Parameters
    ----------
    nvars : int
        The number of variables in the quadratization found
    global_nvars : int
        The minimum number of variables allowed
        
    Returns
    -------
    bool
        True if the number of variables of the quadratization found is greater 
        than the global, False otherwise
    """
    if nvars >= global_nvars: return True 
    return False

def pruning_rule_time(start_time, max_time):
    """Pruning rule based on the time elapsed.
    
    Parameters
    ----------
    start_time : float
        The time when the algorithm started
    max_time : float
        The maximum time allowed
        
    Returns
    -------
    bool
        True if the time elapsed is greater than the maximum time allowed, False otherwise
    """
    if time.time() - start_time > max_time: return True
    return False

def pruning_rule_order(new_vars, max_order):
    """Pruning rule based on the maximum order of derivatives allowed.
    
    Parameters
    ----------
    new_vars : list
        List of proposed new variables
    max_order : int
        The maximum order allowed
        
    Returns
    -------
    bool
        True if the maximum order of the derivatives in the new vars proposed is greater
        than the maximum order allowed, False otherwise
    """
    for var in new_vars:
        if get_diff_order(var) > max_order/2: return True
    return False

def bnb(new_vars, best_nvars, poly_syst, sort_fun, max_der_order):
    """Branch and bound algorithm to find the best quadratization of a polynomial system.
    
    Parameters
    ----------
    new_vars : list
        List of proposed new variables
    best_nvars : int
        The minimum number of variables found in a quadratization
    poly_syst : PolySys
        The polynomial system
    sort_fun : function
        The function to sort the proposed new variables 
        
    Returns
    -------
    tuple
        a tuple with the best quadratization found, the number of variables in the 
        quadratization and the total number of traversed nodes   
    """
    if pruning_rule_nvars(len(new_vars), best_nvars):
        return None, math.inf, 1
    
    if not max_der_order:
        if pruning_rule_order(new_vars, poly_syst.get_max_order()):
                return None, math.inf, 1
    else:
        if pruning_rule_order(new_vars, max_der_order):
                return None, math.inf, 1
    
    poly_syst.set_new_vars(new_vars)
    result_quad = poly_syst.try_make_quadratic()
    
    if result_quad[0]:
        shrinked_quad = shrink_quad(new_vars, poly_syst)
        return shrinked_quad, len(shrinked_quad), 1
    else: 
        if len(new_vars) >= best_nvars - 1:
            return None, math.inf, 1
    
    min_nvars = best_nvars
    best_quad_vars = None
    traversed_total = 1
    prop_vars = prop_new_vars(result_quad[1], new_vars, sort_fun)
    
    for p_vars in prop_vars: 
        quad_vars, nvars, traversed = bnb(new_vars + list(p_vars), min_nvars, poly_syst, sort_fun, max_der_order)
        traversed_total += traversed
        if nvars < min_nvars:
            min_nvars = nvars
            print('Best quadratization until now:', min_nvars, quad_vars)
            best_quad_vars = quad_vars
    
    return best_quad_vars, min_nvars, traversed_total
