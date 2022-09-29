from fractions import Fraction
from functools import reduce
from sympy import *
import numpy as np

def is_a_quadratization(V, deriv):
    V2 = list(set((m1[0]*m2[0], m1[1]*m2[1]) for m1 in V for m2 in V))
    V2_poly = list(term[1] for term in V2)
    quad = []
    
    for pol in deriv:
        if pol[1] not in V2_poly:
            result = is_linear_combination(V2, pol[1])
            if not result: return False
            quad.append(Eq(pol[0], result))
        else: 
            quad.append(Eq(pol[0], V2[V2_poly.index(pol[1])][0]))
    
    print("\nQuadratization:")
    for i in range(len(quad)):
        pprint(exp(quad[i]))
        #quad[i] = poly(quad[i])      
    return quad

def is_linear_combination(V2_names, der_pol):
    V2 = list(map(lambda p: (p[1], p[1].monoms(), p[1].coeffs()), V2_names))
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
            der_expr += s * V2_names[i][0]
            print(f"{s} * {V2_names[i][1]}")
    return simplify(der_expr)
    
# Tests     
u, ux, uxx, uxxx = symbols('u ux uxx uxxx')
w0, w0x, w0xx, w0xxx = symbols('w0 w0x w0xx w0xxx')
u_t, w_0t = symbols('u_t w_0t')

# u_t = u**2 * uxx
# w = u**2
# w_t = 2u * (u**2 * uxx)
V0 = list(map(lambda v, l: (l, poly(v, [u, ux, uxx])), [1, u, ux, uxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx], [1, u, ux, uxx, w0, w0x, w0xx]))
w0t = [(w_0t, poly(2*u**3*uxx, [u, ux, uxx]))]
assert is_a_quadratization(V0, w0t) == [Eq(w_0t, w0*w0xx - w0x**2/2)]

# u_t = u * (ux**2 + u * uxx)
# w = u**2
# w_t = 2u**2 * (ux**2 + u * uxx)
V1 = list(map(lambda v, l: (l, poly(v, [u, ux, uxx])), [1, u, ux, uxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx], [1, u, ux, uxx, w0, w0x, w0xx]))
w0t1 = [(u_t, poly(u*ux**2 + u**2*uxx, [u, ux, uxx])), (w_0t, poly(2*u**2*ux**2 + 2*u**3*uxx, [u, ux, uxx]))]
assert is_a_quadratization(V1, w0t1)#  == [Eq(u_t, u*w0xx/2), Eq(w_0t, w0*w0xx)]

# u_t = u * (2 ux * uxx + u * uxxx + 1)
# w = u**2
# w_t = 2 * u**2 * (2 ux * uxx + u * uxxx + 1)
V2 = list(map(lambda v, l: (l, poly(v, [u, ux, uxx, uxxx])), [1, u, ux, uxx, uxxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx, 4*ux*uxx + 2*u*uxxx],
    [1, u, ux, uxx, uxxx, w0, w0x, w0xx, w0xxx]))
w0t2 = [(u_t, poly(u*(2*ux*uxx + u*uxxx + 1), [u, ux, uxx, uxxx])), (w_0t, poly(2*u**2*(2*ux*uxx + u*uxxx + 1), [u, ux, uxx, uxxx]))]
assert is_a_quadratization(V2, w0t2) #== [Eq(u_t, ux*w0xx/2), Eq(w_0t, w0*w0xx + 2*w0)]
