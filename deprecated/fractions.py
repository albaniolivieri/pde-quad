from sympy import fraction, factor_list, groebner, symbols, simplify, Add, Mul, together
from functools import reduce

def decompose_fraction(pol, pol_vars, i=0):
    """Decomposes a fraction into a polynomial and a list of relations
    
    Parameters
    ----------
    pol : sympy.Expr
        fraction to be decomposed
    pol_vars : list[sympy.PolyElement]
        list of symbols in the expression
    i : int (optional)
        counter for the new symbols
    
    Returns
    -------
    tuple
        a tuple with the new polynomial, the list of relations, the new symbols, 
        the groebner basis and the new counter
    """
    n, d = fraction(together(pol))
    d_factor = factor_list(d)
    
    q_symb, pow_list, fac_list = [], [], []
    for j in range(len(d_factor[1])):
        q_symb.append(symbols(f'q_{i}'))
        pow_list.append(d_factor[1][j][1])
        fac_list.append(d_factor[0] * d_factor[1][j][0])
        i += 1
    
    rel_list = [q * fac - 1 for q, fac in zip(q_symb, fac_list)]
    poly_q = reduce(lambda y, x: y * x[0]**x[1], list(zip(q_symb, pow_list)), 1) * n
    groeb = groebner(rel_list, q_symb + pol_vars, order='lex')
    new_pol = groeb.reduce(poly_q)[1]
    
    return (new_pol, fac_list, q_symb, groeb), i

def get_frac_decomp(exprs, syms):
    """
    Gets the fraction decomposition variables and relations from a PDE
    
    Parameters
    ----------
    exprs : list
        Tuples with the symbol and equations of the PDE
    syms : list
        List of symbols in the PDE
        
    Returns
    -------
    list
        a list with the fraction decompositions
    """
    frac_decomps = []
    i = 0
    flag = False
    for lhs, expr in exprs:
        # Gleb: you kind of assume that you have a sum of fractions as input, is that correct?
        # this may be a bit dangerous
        # `cancel` function from here looks to me like a possible solution https://docs.sympy.org/latest/tutorials/intro-tutorial/simplification.html
        for term in Add.make_args(expr):
            for term2 in Mul.make_args(term):
                if simplify(term2.as_base_exp()[1]).is_negative:
                    flag = True
                    groeb_decomp, i = decompose_fraction(expr, syms, i)
                    # Gleb: this would really benefit from some class holding the data, let us discuss this
                    groeb_decomp = [(lhs, groeb_decomp[0]), groeb_decomp[1], groeb_decomp[2], groeb_decomp[3]]
                    # Gleb: so do you end up having actually several Groebner bases?
                    frac_decomps.append(groeb_decomp)
                    break
            if flag: break
            
    return frac_decomps

# Gleb: do not see this function used
def try_reduce(pol, groeb):
    return groeb.reduce(pol)[1]
    

