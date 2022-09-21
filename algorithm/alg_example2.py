from ast import Eq
from sympy import *

u, ux, uxx = symbols('u ux uxx')

p = Poly(2*ux**2+2*u*uxx)

print(type(p))

V = [1, u, ux, uxx, u**2, u*ux, p]
V2 = []

for var in V: 
    for var2 in V:
        mult = var*var2
        #print(mult)
        if mult not in V2:
            V2.append(mult)
            print(f"{mult} \n")

NS = []

w0t= 2*u**3*uxx

if w0t not in V2: 
    NS.append(w0t)

new_NS = []

for mon in NS: 
    new_NS.append(mon.as_coeff_Mul())
    
print(new_NS)

i = 0
while i < len(V2):
    if type(V2[i]).__name__ == 'Poly':
        #print("aqui", V2[i])
        V2[i] = (V2[i], [prod(x**k for x, k in zip(V2[i].gens, mon)) for mon in V2[i].monoms()], V2[i].coeffs())    
    i += 1  
    
polys = []
for exp in V2: 
    if type(exp) == tuple:
        polys.append(exp)
        print(exp)

#print(polys)

def find_mon(mon, exp_list):
    coinc = []
    for exp in exp_list:
        if type(exp) == tuple:
            if mon in exp[1]: coinc.append(exp)
        else: 
            if mon == exp: coinc.append(exp)
    return coinc

for mon in new_NS:
    occur = find_mon(mon[1], polys)
    index_0 = occur[0][1].index(mon[1])
    occur[0][1].remove(mon[1])

    for mon in occur[0][1]:
        index_1 = occur[0][1].index(mon)
        coinc = find_mon(mon, V2)
        coinc.remove(occur[0])
        for term in coinc:
            poly = occur[0][0]
            if poly - term*occur[0][2][index_1] - NS[0] == 0:
                result = 'quad'
                break
            else: 
                result = 'false'
            
print(result)       
    
#next, how can I avoid an infinite loop 