from functools import reduce
from sympy import *

def is_a_quadratization(V, deriv):
    V2 = list(set((m1[0] * m2[0], m1[1] * m2[1]) for m1 in V for m2 in V))
    V2_poly = [term[1] for term in V2]
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
           
    return quad

def is_linear_combination(V2_names, der_pol):
    V2 = [term[1] for term in V2_names]
    #[print("\nV2 poly", pol) for pol in V2]       
    
    base = list(reduce(lambda base, pol: set(base).union(set(pol.monoms())), V2, []))
    #print(f"\nbase: {base}, length: {len(base)}\n")
    
    lambdas = symbols(["Lambda" + "_%d" % i for i in range(len(V2))])
    subst_lambdas = [(coef, 0) for coef in lambdas]   
    #print(f"lambda set {lambdas}\n")
    
    #print(f"derivative {der_pol.monoms()}\n")
    
    b_vector = zeros(1, len(base), rational=True)
    for monom in der_pol.monoms():
        if monom in base:
            b_vector[base.index(monom)] = Rational(der_pol.coeff_monomial(monom))
        else:
            print("Not a quadratization")
            return False 
    #print(f"b vector {b_vector}\n")
        
    matrix_A = zeros(len(base), len(V2), rational=True)
    for i in range(len(V2)):
        for mon in V2[i].monoms():
            matrix_A[base.index(mon), i] = Rational(V2[i].coeff_monomial(mon))
        
    system = (matrix_A, b_vector)
    #print(f"System: {system}\n")
    sols = list(linsolve(system, lambdas))
    
    if sols == [] or sols[0] == EmptySet:
        print("Not a quadratization")
        return False
    
    sols = list(map(list, sols)) 
    for i in range(len(sols[0])):
        sols[0][i] = sols[0][i].subs(subst_lambdas)

    #print(f"System solution: {sols[0]} \n")
    print("Linear combination:")
            
    der_expr = 0
    for i, s in enumerate(sols[0]):
        if s != 0: 
            der_expr += s * V2_names[i][0]
            print(f"{s} * {V2_names[i][1]}")
            
    return simplify(der_expr)