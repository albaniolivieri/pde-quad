from sympy import *
from sympy import Derivative as D
from .utils import remove_vars   

def prop_new_vars(NS_list, accum_vars, sort_fun):
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
    
    sorted_vars = sort_vars(list_vars, sort_fun)
    return sorted_vars

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

def sort_vars(var_list, fun):
    sorted_list = sorted(var_list, key=fun)
    return sorted_list

def get_order(pol):
    derivs = [x for x in pol.ring.gens if pol.diff(x) != 0]
    order = 0
    for var in derivs:
        if str(var)[-1].isnumeric():  
            order += int(str(var)[-1])         
    return order

def by_degree_order(vars_tup):
    if len(vars_tup) > 1:
        return (max([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]), 
                max([get_order(vars_tup[0]), get_order(vars_tup[1])]))
    else: 
        return (sum(vars_tup[0].degrees()), get_order(vars_tup[0]))
    
def by_order_degree(vars_tup):
    deg, order = 0, 0
    if len(vars_tup) > 1:
        order, deg = (max([get_order(vars_tup[0]), get_order(vars_tup[1])]),
                max([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]))
    else: 
        order, deg = (get_order(vars_tup[0]), sum(vars_tup[0].degrees()))
    return (order, deg)
    
def by_fun(vars_tup):
    if len(vars_tup) > 1:
        return max([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]) + \
                2*max([get_order(vars_tup[0]), get_order(vars_tup[1])])
    else: 
        return sum(vars_tup[0].degrees()) + 2*get_order(vars_tup[0])