from sympy import symbols
from .PolySys import PolySys
from .branch_and_bound import bnb 
from .var_selection import by_fun
from .utils import get_frac_vars

def quadratize(func_eq, n_diff, sort_fun=by_fun, nvars_bound=5):
    """Quadratizes a given PDE
    
    Parameters
    ----------
    func_eq : list[tuple]
        Tuples with the symbol and equations of the PDE
    n_diff : int
        The number of second variable differentiations to do
    sort_fun : function, optional
        The function to sort the proposed new variables
    nvars_bound : int, optional
        The maximum number of variables in the quadratization
    
    Returns
    -------
    tuple
        a tuple with the best quadratization found, the number of variables in the 
        quadratization and the total number of traversed nodes   
    """
    undef_fun = [symbol for symbol, _,  in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    
    _, vars_frac = get_frac_vars(func_eq, undef_fun)
    
    poly_syst = PolySys(func_eq, n_diff, x_var, vars_frac= vars_frac)
    quad = bnb([], nvars_bound, poly_syst, sort_fun)
    
    vars_frac_intro = [1/rel for _, rel in vars_frac]
    
    return quad[0], vars_frac_intro
    

    