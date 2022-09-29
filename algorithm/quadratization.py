from fractions import Fraction
from functools import reduce
from sympy import *
import numpy as np

# maybe V + name_var --> dictionary
def is_a_quadratization(V, deriv, name_var):
    # turn into a dict too
    V2, V2_names = [], []
    for i in range(len(V)):
        for j in range(len(V)):
            mult = V[i]*V[j]
            if mult not in V2:
                V2.append(mult)
                V2_names.append(name_var[i]*name_var[j])
    
    quad = []
    for pol in deriv:
        if pol[1] not in V2:
            result = is_linear_combination(V2, pol[1], V2_names)
            if not result: return False
            quad.append(f"\n{pol[0]} = {result}")
        else: 
            quad.append(f"\n{pol[0]} = {V2_names[V2.index(pol[1])]}")
    
    print("\nQuadratization:")
    [print(new_expr) for new_expr in quad]
    # to return quad but as polynomials
    return True

def is_linear_combination(V2, der_pol, names_V2):
    V2 = list(map(lambda p: (p, p.monoms(), p.coeffs()), V2))
    der_pol = (der_pol, der_pol.monoms(), der_pol.coeffs())
    [print("\nV2 poly", pol) for pol in V2]       
    
    base = list(reduce(lambda base, pol: set(base).union(set(pol[1])), V2, []))
    print(f"\nbase: {base}, length: {len(base)}\n")
    
    lambdas = symbols(["Lambda" + "_%d" % i for i in range(len(V2))])
    subst_lambdas = list((coef, 0) for coef in lambdas)
    
    print(f"lambda set {lambdas}\n")
    
    print(f"derivative {der_pol[1]}\n")
    
    b_vector = zeros(1, len(base), rational=True)
    for i in range(len(der_pol[1])):
        if der_pol[1][i] in base:
            b_vector[base.index(der_pol[1][i])] = Rational(der_pol[2][i])
        else:
            print("Not a quadratization")
            return False 
    print(f"b vector {b_vector}\n")
        
    matrix_A = zeros(len(base), len(V2), rational=True)
    for i in range(len(V2)):
        for j, term in enumerate(V2[i][1]):
            matrix_A[base.index(term), i] = Rational(V2[i][2][j])
        
    system = (matrix_A, b_vector)
    print(f"System: {system}\n")
    sols = list(linsolve(system, lambdas))
    
    if sols == [] or sols[0] == EmptySet:
        print("Not a quadratization")
        return False
    
    sols = list(map(list, sols)) 
    for i in range(len(sols[0])):
        sols[0][i] = sols[0][i].subs(subst_lambdas)

    print(f"System solution: {sols[0]} \n")
    print("Linear combination:")
            
    der_expr = 0
    for i, s in enumerate(sols[0]):
        if s != 0: 
            der_expr += s * names_V2[i]
            print(f"{s} * {V2[i][0]}") 
    print()        
    return simplify(der_expr)
    
# Tests     
u, ux, uxx, uxxx = symbols('u ux uxx uxxx')
w0, w0x, w0xx, w0xxx = symbols('w0 w0x w0xx w0xxx')

# u_t = u**2 * uxx
# w = u**2
# w_t = 2u * (u**2 * uxx)
V0 = list(map(lambda v: poly(v, [u, ux, uxx]), [1, u, ux, uxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx]))
w0t = [("w_t", poly(2*u**3*uxx, [u, ux, uxx]))]
assert is_a_quadratization(V0, w0t, [1, u, ux, uxx, w0, w0x, w0xx]) 

# u_t = u * (ux**2 + u * uxx)
# w = u**2
# w_t = 2u**2 * (ux**2 + u * uxx)
V1 = list(map(lambda v: poly(v, [u, ux, uxx]), [1, u, ux, uxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx]))
w0t1 = [("u_t", poly(u*ux**2 + u**2*uxx, [u, ux, uxx])), ("w_t", poly(2*u**2*ux**2 + 2*u**3*uxx, [u, ux, uxx]))]
assert is_a_quadratization(V1, w0t1, [1, u, ux, uxx, w0, w0x, w0xx]) 

# u_t = u * (2 ux * uxx + u * uxxx + 1)
# w = u**2
# w_t = 2 * u**2 * (2 ux * uxx + u * uxxx + 1)
V2 = list(map(lambda v: poly(v, [u, ux, uxx, uxxx]), [1, u, ux, uxx, uxxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx, 4*ux*uxx + 2*u*uxxx]))
w0t2 = [("u_t", poly(u*(ux**2 + u*uxx + 1), [u, ux, uxx, uxxx])), ("w_t", poly(2*u**2*(2*ux*uxx + u*uxxx + 1), [u, ux, uxx, uxxx]))]
assert is_a_quadratization(V2, w0t2, [1, u, ux, uxx, uxxx, w0, w0x, w0xx, w0xxx]) 
