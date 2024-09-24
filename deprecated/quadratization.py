from sympy import *
from ..algorithm.utils import reduction_sparse

def is_quadratization(V, deriv):
    V2 = list(set((m1[0] * m2[0], m1[1] * m2[1]) for m1 in V for m2 in V))
    V2_poly, names = [], []
    for name, polyn in V2: 
        names.append(name)
        V2_poly.append(polyn)
    quad = []
    V2_red = reduce_set(V2)
    for name, pol in deriv:
        if pol not in V2_poly:
            result = is_linear_combination(V2_red, pol)
            if not result: return False
            else: quad.append(Eq(name, result))
        else: quad.append(Eq(name, names[V2_poly.index(pol)]))
    ("\nQuadratization:")
    for exprs in quad:
        p(exprs)       
    return quad

def reduce_set(V2):
    for i in range(len(V2)):
        for j in range(i):
            V2[i] = reduction_sparse(V2[i], V2[j])
        if V2[i][1] != 0:
            LC = V2[i][1].coeff(V2[i][1].leading_monom())
            V2[i] = (V2[i][0] / LC, V2[i][1] * (1 / LC), V2[i][1].leading_monom())
            for j in range(i):
                V2[j] = reduction_sparse(V2[j], V2[i])                 
    return [(a[0], a[1]) for a in V2]

def is_linear_combination(V2, der_pol):  
    der_tuple = (0, der_pol, der_pol.leading_monom())
    V2 = [(name, pol, pol.leading_monom()) for name, pol in V2]   
    for i in range(len(V2)):
        der_tuple = reduction_sparse(der_tuple, V2[i])
        if der_tuple[1] == 0:
            return simplify(-der_tuple[0])  
    ('rest', der_tuple[1])            
    ("Not a quadratization")
    return False  