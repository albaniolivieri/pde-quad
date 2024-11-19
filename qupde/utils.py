from sympy import symbols, sympify, Eq, Poly, Expr, Float, Rational
from sympy import Derivative as D
from functools import reduce
from itertools import chain, combinations
from fractions import Fraction

def get_order(set_derivs):
    """Gets the maximum order of derivatives in a list of derivatives

    Parameters
    ----------
    set_derivs : list[sympy.Symbol]
        List with the derivatives as expressions

    Returns
    -------
    int
        the maximum order of derivatives
    """
    
    max_order = 0
    for deriv in set_derivs:
        max_order = reduce(max, [der.args[1][1] for der in deriv.atoms(D)], max_order)
    return max_order

def reduction(pol1, pol2):
    """Reduces the first polynomial by the second one

    Parameters
    ----------
    pol1 : sympy.PolyElement
        the polynomial to be reduced
    pol2 : sympy.PolyElement
        the polynomial to reduce by

    Returns
    -------
    sympy.PolyElement
        the reduced polynomial
    """
    coef = pol1[1].coeff_monomial(pol2[2])
    if coef != 0:
        new_pol = pol1[1] - coef * pol2[1]
        return (pol1[0] - coef * pol2[0], new_pol, new_pol.LM())
    else: 
        return pol1

def reduction_sparse(pol1, pol2):
    """Reduces the polynomial pol1 by the polynomial pol2 in the sparse representation

    Parameters
    ----------
    pol1 : tuple
        the polynomial to be reduced as a tuple with the symbol, the polynomial and the leading monomial
    pol2 : tuple
        the polynomial to reduce by as a tuple with the symbol, the polynomial and the leading monomial

    Returns
    -------
    tuple
        the reduced polynomial as a tuple with the symbol, the polynomial and the leading monomial
    """
    if pol2[1] != 0:
        coef = pol1[1].coeff(pol2[2])
        if coef != 0:
            new_pol = pol1[1] - coef * pol2[1]
            return (pol1[0] - coef.as_expr() * pol2[0], new_pol, new_pol.leading_monom()) 
    return pol1

def diff_dict(pol, dic, frac_decomp=None, order=1):
    """differentiates a polynomial by a dictionary of variables

    Parameters
    ----------
    pol : sympy.PolyElement
        polynomial to be differentiated
    dic : dict
        dictionary with the variables to differentiate by
    order : int, optional
        the order of differentiation (default is 1)

    Returns
    -------
    sympy.PolyElement
        the differentiated polynomial
    """
    deriv = pol
    for _ in range(1, order + 1):
        deriv = sum(deriv.diff(k) * v for (k, v) in dic.items())
    if frac_decomp:
        return pol.ring(frac_decomp.try_reduce(deriv))
    return deriv

def remove_vars(list_vars, accum_vars, axis):
    """Removes variables already proposed in the quadratization or of degree less than two

    Parameters
    ----------
    list_vars : list[sympy.PolyRing]
        Set of variables to check
    accum_vars : list[sympy.PolyRing]
        List with all variables proposed up to this point
    axis : int
        Axis to check
    Returns
    -------
    list[sympy.PolyRing]
        the filtered list of variables
    """
    for i in range(len(list_vars)):
        if len(list_vars[i]) > 1:
            if list_vars[i][axis] in accum_vars or (sum(list_vars[i][axis].degrees()) <= 1):
                list_vars[i] = (list_vars[i][int(not axis)],)
            elif list_vars[i][axis] == list_vars[i][int(not axis)]:
                list_vars[i] = (list_vars[i][axis],)
    return list_vars

def powerset(iterable):
    """Returns the powerset of a set

    Parameters
    ----------
    iterable : list[sympy.PolyRing]
        Set to get the powerset from

    Returns
    -------
    list[sympy.PolyRing]
        the powerset of the set
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)))

def shrink_quad(quad_vars, poly_syst):
    """Checks if the quadratization can be shrunk to a smaller set of variables.
    
    Parameters
    ----------
    quad_vars : list
        List of variables in the quadratization
    poly_syst : RatSys
        The polynomial system
        
    Returns
    -------
    list
        a list with a quadratization of an equal or lesser order than the original
    """
    final_vars = quad_vars
    subsets = powerset(quad_vars)
    for var_group in subsets: 
        poly_syst.set_new_vars(var_group)
        res, _ = poly_syst.try_make_quadratic() 
        if res:
            return list(var_group)
    return final_vars

def ring_to_expr(ring_pol, ring_syms=[], consts=[]):
    """converts a polynomial ring to a sympy expression

    Parameters
    ----------
    ring_pol : sympy.PolyRing
        polynomial to be converted
    ring_syms : list[sympy.PolyElement]
        symbols of the polynomial ring
    consts : list[sympy.PolyElement]
        constants of the polynomial ring

    Returns
    -------
    tuple
        a tuple with a the polynomial as a sympy expression and the symbols as a list of sympy symbols
    """
    expr_syms = {}
    str_consts = [const.name for const in consts]
    if ring_syms:
        for var in ring_syms:
            if str(var) in str_consts:
                expr_syms[str(var)] = consts[str_consts.index(str(var))]
            else: 
                expr_syms[str(var)] = symbols(str(var))
        expr_pol = sympify(str(ring_pol), locals=expr_syms) 
    else:
        expr_pol = sympify(str(ring_pol))
    return expr_pol

def expr_to_ring(R, expr_pol):
    """converts a sympy expression to a polynomial ring

    Parameters
    ----------
    R : sympy.PolyRing
        Polynomial ring to convert to
    expr_pol : sympy.Expr
        expression to be converted

    Returns
    -------
    sympy.PolyRing
        the resulting polynomial ring
    """
    return R(expr_pol)         

def get_diff_order(pol):
    """Returns the order of the highest derivative in a polynomial
    
    Parameters
    ----------
    pol : sympy.PolyElement
        Polynomial to check
    
    Returns
    -------
    int
        the order of the highest derivative in the polynomial
    """
    derivs = [x for x in pol.ring.gens if pol.diff(x) != 0]
    order = 0
    for var in derivs:
        if str(var)[-1].isnumeric():  
            order += int(str(var)[-1])         
    return order
