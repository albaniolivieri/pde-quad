from functools import reduce
import sympy as sym

# from sympy import symbols, Add, Mul, simplify
from sympy import Derivative as D
from .PolySys import *
from .fractions import decompose_fraction


def test_quadratization(func_eq, new_vars: list, n_diff: int):
    undef_fun = [symbol for symbol, _ in func_eq]
    x_var = [
        symbol for symbol in undef_fun[0].free_symbols if symbol != sym.symbols("t")
    ].pop()

    print("func_eq", func_eq)

    symbols = [undef_fun]
    frac_rel = []
    vars_frac = []
    i = 0

    for _, expr in func_eq:
        print("expression", expr)
        for term in sym.Add.make_args(expr):
            for term2 in sym.Mul.make_args(term):
                print("exponent", term2.as_base_exp()[1])
                if sym.simplify(term2.as_base_exp()[1]).is_negative:
                    frac_decomp = decompose_fraction(term2, undef_fun)
                    print("fraction decomposition", frac_decomp)
                    for rel in frac_decomp[1]:
                        if rel not in frac_rel: 
                            frac_rel.append(rel)
                            vars_frac.append((sym.symbols(f'q{i}'), rel))
                            i += 1
                    
    print("vars_frac", vars_frac)

    # TODO: check if there's a fraction in func_eq => send a flag

    poly_syst = PolySys(func_eq, n_diff, x_var, new_vars, vars_frac)
    
    # convert back before returning
    # return poly_syst.try_make_quadratic()