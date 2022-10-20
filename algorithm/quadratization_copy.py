from functools import reduce
import time
from sympy import *

def is_a_quadratization(V, deriv):
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
            else: 
                quad.append(Eq(name, result))
        else: 
            quad.append(Eq(name, names[V2_poly.index(pol)]))
    
    print("\nQuadratization:")
    for exprs in quad:
        pprint(exprs)  
           
    return quad

def reduce_set(V2):
    print(len(V2))
    ti = time.time()
    for i in range(len(V2)):
        for j in range(i):
            coef = V2[i][1].coeff_monomial(V2[j][2])
            if coef != 0:
                V2[i] = (V2[i][0] - coef * V2[j][0], V2[i][1] - coef * V2[j][1])
        LC = V2[i][1].LC()
        V2[i] = (V2[i][0] / LC, V2[i][1] * (1 / LC), V2[i][1].LM())
        for j in range(i):
            coef = V2[j][1].coeff_monomial(V2[i][2])
            if coef != 0:
               V2[j] = (V2[j][0] - coef * V2[i][0], V2[j][1] - coef * V2[i][1], V2[j][2])
    print(time.time()-ti)                     
    return [(a[0], a[1]) for a in V2]

def is_linear_combination(V2, der_pol):          
    der_expr = 0
    for name_1, pol1 in V2:
        subt = der_pol
        if der_pol.LM() in pol1.monoms():
            LM_subt = subt.LM()
            LC_subt = subt.LC()
            coef = pol1.coeff_monomial(LM_subt)
            subt -= (LC_subt/coef)*pol1
            der_expr += (LC_subt/coef)*name_1
            if subt == 0:
                return simplify(der_expr)
            for name, pol2 in V2:
                LM_subt = subt.LM()
                LC_subt = subt.LC()
                if LM_subt in pol2.monoms():
                    coef = pol2.coeff_monomial(LM_subt)
                    subt -= (LC_subt/coef)*pol2
                    der_expr += (LC_subt/coef)*name
                    if subt == 0:
                        return simplify(der_expr)                 
    print("Not a quadratization")
    return False            

