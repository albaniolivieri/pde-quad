from functools import reduce
from sympy import *

def is_a_quadratization(V, deriv):
    V2 = list(set((m1[0] * m2[0], m1[1] * m2[1]) for m1 in V for m2 in V))
    V2_poly, names = [], []
    for name, polyn in V2: 
        names.append(name)
        V2_poly.append(polyn)
    quad = []
    
    base, matrix_A, lambdas = get_matrix_system(V2)
    
    for pol in deriv:
        if pol[1] not in V2_poly:
            result = is_linear_combination(names, base, matrix_A, lambdas, pol[1])
            if not result: return False
            quad.append(Eq(pol[0], result))
        else: 
            quad.append(Eq(pol[0], V2[V2_poly.index(pol[1])][0]))
    
    print("\nQuadratization:")
    for exprs in quad:
        pprint(exprs)  
           
    return quad

def is_linear_combination(names, base, matrix_A, lambdas, der_pol):
    b_vector = zeros(len(base), 1, rational=True)
    for monom in der_pol.monoms():
        if monom in base:
            b_vector[base.index(monom)] = Rational(der_pol.coeff_monomial(monom))
        else:
            print("Not a quadratization")
            return False 
        
    system = (matrix_A, b_vector)
    sols = list(linsolve(system, lambdas[0]))
    
    if sols == [] or sols[0] == EmptySet:
        print("Not a quadratization")
        return False
    
    sols = list(map(list, sols)) 
    for i in range(len(sols[0])):
        sols[0][i] = sols[0][i].subs(lambdas[1])
            
    der_expr = 0
    for i, s in enumerate(sols[0]):
        if s != 0: 
            der_expr += s * names[i]
            
    return simplify(der_expr)

def get_matrix_system(V2_names):
    V2 = [term[1] for term in V2_names]
    base = list(reduce(lambda base, pol: set(base).union(set(pol.monoms())), V2, []))
    
    matrix_A = zeros(len(base), len(V2), rational=True)
    for i in range(len(V2)):
        for mon in V2[i].monoms():
            matrix_A[base.index(mon), i] = Rational(V2[i].coeff_monomial(mon))
            
    lambdas = symbols(["Lambda" + "_%d" % i for i in range(len(V2))])
    subst_lambdas = [(coef, 0) for coef in lambdas] 
    
    return base, matrix_A, (lambdas, subst_lambdas)
