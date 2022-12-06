from sympy import *
from sympy import Derivative as D

def propose_variables(func_eq):
    undef_fun = [symbol for symbol, _ in func_eq]
    

def sort_NS(ns_pol):
    for monom in ns_pol.itermonoms():
        list_vars = divide_mon(ns_pol, monom)
        by_degree_order = sort_degree_order(list_vars)
    return by_degree_order
        
def divide_mon(ns_pol, mon):
    prop_vars = []
    prop_vars.append(ns_pol.ring({mon:1}))
    if len(set(mon)) == 1 and mon[0]%2 == 0:
        prop_vars.append(ns_pol.ring({tuple(int(item/2) for item in mon):1}))
    for i in range(len(mon)):
        aux = list(mon)
        if aux[i] != 0:
            aux[i] = aux[i] - 1
            prop_vars.append(ns_pol.ring({tuple(aux):1}))
    print('prop_vars', prop_vars)
    return prop_vars

def get_order(pol):
    derivs = [x for x in pol.ring.gens if pol.diff(x) != 0]
    order = 0
    for var in derivs:
        if str(var)[-1].isnumeric():  
            order = max(int(str(var)[-1]), order)                
    return order

def sort_degree_order(var_list):
    sort_degree_order = sorted(var_list, key=lambda x: (sum(x.degrees()), get_order(x)))
    print('sort_degree', sort_degree_order)
    return sort_degree_order

R, u, u_x1, u_x2, u_x3 = ring('u, u_x1, u_x2, u_x3', QQ)
#ut5 = u**2*ux + u
#print(selection(ut5))

ut5 = u**3 * u_x3
print(selection(ut5))