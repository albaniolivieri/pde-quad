from sympy import *
from sympy import Derivative as D
from functools import reduce
from itertools import chain, combinations

def get_order(set_derivs):
    max_order = 0
    for deriv in set_derivs:
        max_order = reduce(max, [der.args[1][1] for der in deriv.atoms(D)], max_order)
    return max_order

def reduction(pol1, pol2):
    coef = pol1[1].coeff_monomial(pol2[2])
    if coef != 0:
        new_pol = pol1[1] - coef * pol2[1]
        return (pol1[0] - coef * pol2[0], new_pol, new_pol.LM())
    else: 
        return pol1

# Gleb: quick docstrings here would be great as well
def reduction_sparse(pol1, pol2):
    if pol2[1] != 0:
        coef = pol1[1].coeff(pol2[2])
        if coef != 0:
            new_pol = pol1[1] - coef * pol2[1]
            return (pol1[0] - coef * pol2[0], new_pol, new_pol.leading_monom())
    return pol1

def diff_dict(pol, dic, order=1):
    deriv = pol
    for _ in range(1, order + 1):
        deriv = sum(deriv.diff(k) * v for (k, v) in dic.items())
    return deriv

def remove_vars(list_vars, accum_vars, axis):
    for i in range(len(list_vars)):
        if len(list_vars[i]) > 1:
            if list_vars[i][axis] in accum_vars or (sum(list_vars[i][axis].degrees()) <= 1):
                list_vars[i] = (list_vars[i][int(not axis)],)
    return list_vars

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)))

def ring_to_expr(ring_syms, ring_pol):
    expr_syms = [symbols(str(var)) for var in ring_syms]
    expr_pol = sympify(str(ring_pol)) 
    return expr_pol, expr_syms

def expr_to_ring(R, expr_pol):
    return R.ring_new(expr_pol)         
