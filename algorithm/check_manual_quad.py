from sympy import symbols

from .PolySys import PolySys

def test_quadratization(func_eq, new_vars: list, n_diff: int, first_indep=symbols('t'), frac_vars=[]):
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
    undef_fun = [symbol for symbol, _,  in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != first_indep].pop()
    frac_vars = [expr for _, expr in frac_vars]

    poly_syst = PolySys(func_eq, n_diff, (first_indep, x_var), new_vars, frac_vars)    
    
    return poly_syst.try_make_quadratic()