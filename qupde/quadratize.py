from sympy import symbols
from .RatSys import RatSys
from .branch_and_bound import bnb
from .nearest_neighbor import nearest_neighbor
from .var_selection import by_fun


def quadratize(
    func_eq,
    n_diff,
    sort_fun=by_fun,
    nvars_bound=10,
    first_indep=symbols("t"),
    max_der_order=None,
    search_alg = 'bnb' # 'bnb' or 'nn'
):
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
    first_indep : sympy.Symbol, optional
        The first independent variable of the PDE
    max_der_order : int, optional
        The maximum order of derivatives allowed in the new variables

    Returns
    -------
    tuple
        a tuple with the best quadratization found, the variables introduced
        from rational functions and the total number of traversed nodes
    """
    undef_fun = [symbol for symbol, _, in func_eq]
    x_var = [
        symbol for symbol in undef_fun[0].free_symbols if symbol != first_indep
    ].pop()

    poly_syst = RatSys(func_eq, n_diff, (first_indep, x_var))
    vars_frac_intro = poly_syst.get_frac_vars()
    
    if search_alg == 'nn':
        quad = nearest_neighbor(poly_syst, sort_fun, new_vars=[])
        if not quad[0] and not vars_frac_intro:
            print("Quadratization not found")
        return quad[0], vars_frac_intro, quad[1]
    elif search_alg == 'bnb':
        quad = bnb([], nvars_bound, poly_syst, sort_fun, max_der_order)
        if not quad[0] and not vars_frac_intro:
            print("Quadratization not found")
        return quad[0], vars_frac_intro, quad[2]

