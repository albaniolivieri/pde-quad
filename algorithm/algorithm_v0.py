from ast import Eq
from sympy import *

def build_V2(V):
    V2 = []
    for var in V: 
        for var2 in V:
            mult = var*var2
            if mult not in V2: V2.append(mult)
    return V2

def build_NS(ders_t, V2):
    NS = []
    for deriv in ders_t:
        if deriv not in V2: 
            NS.append((deriv, deriv.as_coeff_Mul()))

def modify_poly(V2):
    i = 0
    while i < len(V2):
        if type(V2[i]).__name__ == 'Poly':
            V2[i] = (V2[i], [prod(x**k for x, k in zip(V2[i].gens, mon)) for mon in V2[i].monoms()], V2[i].coeffs())    
        i += 1  

def find_mon(mon, exp_list):
    coinc = []
    for exp in exp_list:
        if type(exp) == tuple:
            if mon in exp[1]: coinc.append(exp)
        else: 
            if mon == exp: coinc.append(exp)
    return coinc

def is_quadratized(V2, NS):
    polys = []
    for exp in V2: 
        if type(exp) == tuple:
            polys.append(exp)
            print(exp)
            
    for mon in NS:
        occur = find_mon(mon[1][1], polys)
        #index_0 = occur[0][1].index(mon[1][1])
        occur[0][1].remove(mon[1][1])

        for mon_p in occur[0][1]:
            index_coef = occur[0][1].index(mon_p)
            coinc = find_mon(mon_p, V2)
            coinc.remove(occur[0])
            for term in coinc:
                poly = occur[0][0]
                if poly - term*occur[0][2][index_coef] - mon_p[0] == 0:
                    result = 'quad'
                    break
                else: 
                    result = 'false'
            
       
    
#next, how can I avoid an infinite loop 