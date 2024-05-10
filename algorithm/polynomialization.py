from sympy import symbols, cancel, exp, sin, cos, Pow

non_pol_funcs = [exp, sin, cos, Pow]

def polynomialize(pde_sys, first_indep=symbols('t')):
    pde_sys = [(lhs, cancel(rhs)) for lhs, rhs in pde_sys]
    
def replace_non_pol_funcs(expr):
    pass

    
    