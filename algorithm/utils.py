from sympy import symbols, sympify, Eq, Add, Mul, simplify
from sympy import Derivative as D
from functools import reduce
from itertools import chain, combinations
from .fractions import decompose_fraction

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
    """Reduces the first polynomial by the second one using sparse polynomials

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
            return (pol1[0] - coef * pol2[0], new_pol, new_pol.leading_monom())
    return pol1

def diff_dict(pol, dic, order=1):
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

def ring_to_expr(ring_syms, ring_pol):
    """converts a polynomial ring to a sympy expression

    Parameters
    ----------
    ring_syms : list[sympy.PolyElement]
        symbols of the polynomial ring
    ring_pol : sympy.PolyRing
        polynomial to be converted

    Returns
    -------
    tuple
        a tuple with a the polynomial as a sympy expression and the symbols as a list of sympy symbols
    """
    expr_syms = [symbols(str(var)) for var in ring_syms]
    expr_pol = sympify(str(ring_pol)) 
    return expr_pol, expr_syms

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
    return R.ring_new(expr_pol)         

def revert_frac_decomp(quad_exprs, frac_rel, symbols_in_ring):
    """Reverts a quadratization without the fraction decomposition variables

    Parameters
    ----------
    quad_exprs : list[sympy.Eq]
        List with the expressions of the quadratization
    frac_rel : list[tuple]
        List of tuples with the relations of the fraction decomposition

    Returns
    -------
    list[sympy.Eq]
        the original quadratization without the fraction decomposition variables
    """
    frac_subs = [(var, 1/ring_to_expr(symbols_in_ring, rel)[0]) for var, rel in frac_rel]
    for i in range(len(quad_exprs)):
        quad_exprs[i] = Eq(quad_exprs[i].lhs, quad_exprs[i].rhs.subs(frac_subs))
    return quad_exprs


def diff_frac(num, den, symb_var, dic):
    """Differentiates a fraction by a dictionary of variables

    Parameters
    ----------
    num : sympy.PolyRing
        numerator of the fraction
    den : sympy.PolyRing
        denominator of the fraction
    symb_var : sympy.PolyRing
        polinomial variable of the fraction
    dic : dict
        dictionary with the variables to differentiate by

    Returns
    -------
    sympy.PolyRing
        the differentiated fraction
    """
    return (diff_dict(num, dic) * den - diff_dict(den, dic) * num) * symb_var**2 

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

def get_frac_vars(func_eq, undef_fun):
    """
    Gets the fraction decomposition variables and relations from a PDE
    
    Parameters
    ----------
    func_eq : list[tuple]
        Tuples with the symbol and equations of the PDE
        undef_fun : list
        List of undefined functions in the PDE
        
    Returns
    -------
    tuple
        tuple with two lists: one with the fraction decomposition relations and
        the other with the fraction decomposition variables
    """
    frac_rel, vars_frac = [], []
    i = 0
    for _, expr in func_eq:
        for term in Add.make_args(expr):
            for term2 in Mul.make_args(term):
                if simplify(term2.as_base_exp()[1]).is_negative:
                    frac_decomp = decompose_fraction(term2, undef_fun)
                    for rel in frac_decomp[1]:
                        if rel not in frac_rel: 
                            frac_rel.append(rel)
                            vars_frac.append((symbols(f'q{i}'), rel))
                            i += 1
    return frac_rel, vars_frac
