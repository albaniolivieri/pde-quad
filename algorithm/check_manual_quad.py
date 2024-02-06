from sympy import symbols

from .PolySys import PolySys
from .utils import get_frac_vars

def test_quadratization(func_eq, new_vars: list, n_diff: int):
    """Test the proposed quadratization of a given PDE
    
    Parameters
    ----------
    func_eq : list[tuple]
        Tuples with the symbol and equations of the PDE
    new_vars : list
        List of proposed new variables
    n_diff : int
        The number of second variable differentiations to do
        
    Returns
    -------
    bool
        True if the proposed quadratization is valid, False otherwise
    """
    undef_fun = [symbol for symbol, _ in func_eq]
    x_var = [
        symbol for symbol in undef_fun[0].free_symbols if symbol != symbols("t")
    ].pop()

    frac_rel, vars_frac = get_frac_vars(func_eq, undef_fun)
    
    if len(frac_rel) > 0:
        frac_subs = [(1/rel, var) for var, rel in vars_frac]
        new_vars = [var.subs(frac_subs) for var in new_vars] 

    poly_syst = PolySys(func_eq, n_diff, x_var, new_vars, vars_frac)    
    
    return poly_syst.try_make_quadratic()