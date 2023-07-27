from sympy import symbols
from .PolySys import PolySys
from .branch_and_bound import bnb 
from .var_selection import by_fun

def quadratize(func_eq, n_diff, sort_fun=by_fun, nvars_bound=10):
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    
    poly_syst = PolySys(func_eq, n_diff, x_var)
    quad = bnb([], nvars_bound, poly_syst, sort_fun)
    
    return quad
    

    