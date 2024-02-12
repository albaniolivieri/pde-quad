from sympy import fraction, factor_list, groebner, symbols
from functools import reduce

def decompose_fraction(frac, pol_vars):
    """Decomposes a fraction into a polynomial and a list of factors
    
    Parameters
    ----------
    frac : sympy.Expr
        fraction to be decomposed
    pol_vars : list[sympy.PolyElement]
        list of symbols in the polynomial ring
    
    Returns
    -------
    tuple
        a tuple with the new polynomial, the list of relations and the new symbols
    """
    n, d = fraction(frac)
    d_factor = factor_list(d)
    
    q_symb, pow_list, fac_list = [], [], []
    for i in range(len(d_factor[1])):
        q_symb.append(symbols(f'q_{i}'))
        pow_list.append(d_factor[1][i][1])
        fac_list.append(d_factor[1][i][0])
    
    rel_list = [q*fac - 1 for q, fac in zip(q_symb, fac_list)]
    poly_q = reduce(lambda y, x: y*x[0]**x[1], list(zip(q_symb, pow_list)), 1) * n
    groeb = groebner(rel_list, q_symb + pol_vars,  order='lex')
    new_pol = groeb.reduce(poly_q)[1]
    
    return new_pol, fac_list, q_symb

