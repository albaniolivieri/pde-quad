from fractions import Fraction
from sympy import *
import numpy as np

def is_a_quadratization(V, deriv):
    V2 = list(set([(m1 * m2) for m1 in V for m2 in V]))
    for mon in deriv:
        if mon not in V2: 
            if not is_linear_combination(V2, mon): return False
    print("\nIt is a quadratization")
    return True

def is_linear_combination(V2, mon):
    V2 = list(map(lambda p: (p, p.monoms(), p.coeffs()), V2))
    mon = (mon, mon.monoms(), mon.coeffs())
    [print("V2 poly", pol) for pol in V2]   
    
    base = []
    for pol in V2:
        for terms_pol in pol[1]:
            if terms_pol not in base:
                base.append(terms_pol)
    print(f"\nbase {base}\n")
    
    lambdas = symbols(["Lambda" + "_%d" % i for i in range(len(V2))])
    print(f"lambda set {lambdas}\n")
    
    b_vector = np.zeros(len(base), dtype=int) + Fraction()
    print(f"derivative mon {mon[1][0]}\n")
    if mon[1][0] in base: 
        b_vector[base.index(mon[1][0])] = Fraction(mon[2][0])
        print(f"b vector {b_vector}\n")
    else: 
        print("Not a quadratization")
        return False
        
    matrix_A = []
    for pol in V2:
        sub_v = np.zeros(len(base), dtype=int) + Fraction()
        for term in pol[1]:
            if term in base: sub_v[base.index(term)] = Fraction(pol[2][pol[1].index(term)])
        matrix_A.append(sub_v)
        
    system = (Matrix(matrix_A).T, Matrix(b_vector))
    print(f"System: {system}\n")
    sols = list(linsolve(system, lambdas))
    
    if sols == EmptySet:
        print("Not a quadratization")
        return False
    
    print(f"System solution: {sols} \n")
    print("Linear combination:")
    i = 0
    while i < len(sols[0]):
        if sols[0][i] != 0: print(f"{sols[0][i]} * {V2[i][0]}") 
        i += 1
    return True
    
# Tests     
u, ux, uxx = symbols('u ux uxx')

V = list(map(lambda v: poly(v, [u, ux, uxx]), [1, u, ux, uxx, u**2, u*ux, 2*ux**2+2*u*uxx]))
w0t= [poly(2*u**3*uxx, [u, ux, uxx])]

V1 = list(map(lambda v: poly(v, [u, ux]), [1, u, ux, u + ux**2]))
w0t1 = [poly(u*ux**2, [u, ux])]

V2 = list(map(lambda v: poly(v, [u, ux]), [1, u, ux, u**2, 2*u*ux]))
w0t2 = [poly(2*u**3*ux, [u, ux]), poly(2*u**4, [u, ux])]

assert is_a_quadratization(V, w0t) 
assert is_a_quadratization(V1, w0t1)
assert is_a_quadratization(V2, w0t2)  



    