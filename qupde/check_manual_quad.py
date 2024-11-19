from sympy import symbols

from .RatSys import RatSys

def check_quadratization(func_eq, new_vars: list, n_diff: int, first_indep=symbols('t')):
    """Checks if a given set of new variables is a quadratization for the provided PDE
    
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
    undef_fun = [symbol for symbol, _,  in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != first_indep].pop()

    poly_syst = RatSys(func_eq, n_diff, (first_indep, x_var), new_vars)    
    
    return poly_syst.try_make_quadratic()