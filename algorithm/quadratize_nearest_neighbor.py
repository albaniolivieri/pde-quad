from sympy import symbols
from .PolySys import PolySys
from .var_selection import by_fun
from .nearest_neighbor import nearest_neighbor

# Gleb: the code seems to be almost the same as in `quadratize`, so maybe it should be a function
# taking a method (bnb or nn) as a parameter?
def quadratize(func_eq, n_diff, sort_fun=by_fun, first_indep=symbols('t')):
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
        a tuple with the best quadratization found, the variables introduced 
        from rational functions and the total number of traversed nodes   
    """
    undef_fun = [symbol for symbol, _, in func_eq]
    x_var = [
        symbol for symbol in undef_fun[0].free_symbols if symbol != first_indep].pop()

    poly_syst = PolySys(func_eq, n_diff, (first_indep, x_var))
    quad = nearest_neighbor(poly_syst, sort_fun, new_vars=[])
    vars_frac_intro = poly_syst.get_frac_vars()

    if not quad[0] and not vars_frac_intro:
        print("Quadratization not found")
    return quad[0], vars_frac_intro, quad[1]
