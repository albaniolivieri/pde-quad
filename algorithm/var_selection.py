from sympy import *
from sympy import Derivative as D
from .utils import remove_vars

def propose_variables(NS_list):
    vars_by_pol = [sort_new_vars(pol[1]) for pol in NS_list]
    prop_vars = list(map(list, zip(*vars_by_pol))) # different sizes in prop vars list eliminates last candidate 
    return prop_vars    

def prop_new_vars(NS_list, accum_vars):
    list_vars = []
    for ns_pol in NS_list:
        for monom in ns_pol[1].itermonoms(): 
            list_vars += get_decompositions(monom)      
    list_vars = list(map(lambda x: (NS_list[0][1].ring({x[0]:1}), NS_list[0][1].ring({x[1]:1})), list_vars))
    
    list_vars = remove_vars(list_vars, accum_vars, 0)
    list_vars = remove_vars(list_vars, accum_vars, 1)
    for i in range(len(list_vars)): 
        if not list_vars[i]: 
            list_vars.remove(list_vars[i])
    
    by_degree_order = sort_order_degree_fun(list_vars)
    print('by_degree_order', by_degree_order)
    return by_degree_order

def get_decompositions(monomial):
    if len(monomial) == 0:
        return {(tuple(), tuple())}
    result = set()
    prev_result = get_decompositions(tuple(monomial[:-1]))
    for r in prev_result:
        for i in range(monomial[-1] + 1):
            a, b = tuple(list(r[0]) + [i]), tuple(list(r[1]) + [monomial[-1] - i])
            result.add((min(a, b), max(a, b)))
    return result

def divide_mon(ns_pol, mon):
    prop_vars = []
    prop_vars.append(ns_pol.ring({mon:1}))
    if (len(set(mon)) == 1 or (len(set(mon)) == 2 and 0 in set(mon))) and mon[0]%2 == 0:
        prop_vars.append(ns_pol.ring({tuple(int(item/2) for item in mon):1}))
    for i in range(len(mon)):
        aux = list(mon)
        if aux[i] != 0:
            aux[i] = aux[i] - 1
            prop_vars.append(ns_pol.ring({tuple(aux):1}))
    return prop_vars

def get_order(pol):
    derivs = [x for x in pol.ring.gens if pol.diff(x) != 0]
    order = 0
    for var in derivs:
        if str(var)[-1].isnumeric():  
            order += int(str(var)[-1])         
    return order

def sorting_fun(vars_tup):
    if len(vars_tup) > 1:
        return (sum([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]), 
                sum([get_order(vars_tup[0]), get_order(vars_tup[1])]))
    else: 
        return (sum(vars_tup[0].degrees()), get_order(vars_tup[0]))
    
def sorting_fun_inv(vars_tup):
    deg, order = 0, 0
    if len(vars_tup) > 1:
        order, deg = (max([get_order(vars_tup[0]), get_order(vars_tup[1])]),
                max([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]))
    else: 
        order, deg = (get_order(vars_tup[0]), sum(vars_tup[0].degrees()))
    return (order, deg)
    
def sorting_degree_order(vars_tup):
    if len(vars_tup) > 1:
        return sum([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]) + \
                2*sum([get_order(vars_tup[0]), get_order(vars_tup[1])])
    else: 
        return sum(vars_tup[0].degrees()) + 2*get_order(vars_tup[0])

def sort_degree_order(var_list):
    sort_degree_order = sorted(var_list, key=sorting_fun)
    return sort_degree_order

def sort_order_degree(var_list):
    sort_order_degree = sorted(var_list, key=sorting_fun_inv)
    return sort_order_degree

def sort_order_degree_fun(var_list):
    sort_order_degree = sorted(var_list, key=sorting_degree_order)
    return sort_order_degree
