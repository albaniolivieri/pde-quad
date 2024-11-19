from .utils import remove_vars, get_diff_order

def prop_new_vars(NS_list, accum_vars, sort_fun):
    """Proposes new variables for the quadratization of a PDE
    
    Parameters
    ----------
    NS_list : list[tuple]
        List of tuples with the symbol and expressions not quadratized of the PDE
    accum_vars : list
        List with all variables proposed up to this point
    sort_fun : function
        Function to sort the variables
        
    Returns
    -------
    list
        the proposed new variables
    """
    list_vars = []
    NS_list = sorted(NS_list, key=lambda x: str(x[0]))
    for ns_pol in NS_list:
        for monom in ns_pol[1].itermonoms(): 
            list_vars += get_decompositions(monom)
            break
        break
                
    list_vars = list(map(lambda x: (NS_list[0][1].ring({x[0]:1}), NS_list[0][1].ring({x[1]:1})), list_vars))
    
    list_vars = remove_vars(list_vars, accum_vars, 0)
    list_vars = remove_vars(list_vars, accum_vars, 1)

    for i in range(len(list_vars)): 
        if not list_vars[i]: 
            list_vars.remove(list_vars[i])
    
    sorted_vars = sort_vars(list_vars, sort_fun)
    return sorted_vars

def get_decompositions(monomial):
    """ Returns the decompositions of a monomial to turn it into linear or quadratic
    
    Parameters
    ----------
    monomial : tuple
        Monomial to decompose
    
    Returns
    -------
    set
        the decompositions of the monomial
    """
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
    """Sorts the list of variables according to the function fun
    
    Parameters
    ----------
    var_list : list
        List of variables to sort
    fun : function
        Function to sort the variables
        
    Returns
    -------
    list
        the sorted list of variables
    """
    sorted_list = sorted(var_list, key=fun)
    return sorted_list

def by_degree_order(vars_tup):
    """Function to sort the variables by the sum of their degrees and the 
    maximum order of differentiation
    
    Parameters
    ----------
    vars_tup : tuple
        Tuple with the variables to be sorted
    
    Returns
    -------
    tuple
        tuple with the sorting criteria given by the sum of the degrees and 
        the maximum order of differentiation of the variables
    """
    if len(vars_tup) > 1:
        return (max([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]), 
                max([get_diff_order(vars_tup[0]), get_diff_order(vars_tup[1])]))
    else: 
        return (sum(vars_tup[0].degrees()), get_diff_order(vars_tup[0]))
    
def by_order_degree(vars_tup):
    """Function to sort the variables by the maximum order of differentiation and the
    sum of their degrees
    
    Parameters
    ----------
    vars_tup : tuple
        Tuple with the variables to be sorted
        
    Returns
    -------
    tuple
        tuple with sorting criteria given by the maximum order of differentiation and 
        the sum of the degrees of the variables
    """
    deg, order = 0, 0
    if len(vars_tup) > 1:
        order, deg = (max([get_diff_order(vars_tup[0]), get_diff_order(vars_tup[1])]),
                max([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]))
    else: 
        order, deg = (get_diff_order(vars_tup[0]), sum(vars_tup[0].degrees()))
    return (order, deg)
    
def by_fun(vars_tup):
    """Function to sort the variables by the function: degree + 2 * order
    
    Parameters
    ----------
    vars_tup : tuple
        Tuple with the variables to be sorted
    
    Returns
    -------
    int
        sorting criteria given by the value of the function degree + 2 * order
    """
    if len(vars_tup) > 1:
        return max([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]) + \
                2*max([get_diff_order(vars_tup[0]), get_diff_order(vars_tup[1])])
    else: 
        return sum(vars_tup[0].degrees()) + 2*get_diff_order(vars_tup[0])
    
def by_fun2(vars_tup):
    """
    Function to sort the variables by the function: degree + 4 * order
    
    Parameters
    ----------
    vars_tup : tuple
        Tuple with the variables to be sorted
    
    Returns
    -------
    int
        sorting criteria given by the value of the function degree + 4 * order
    """
    if len(vars_tup) > 1:
        return sum([sum(vars_tup[0].degrees()), sum(vars_tup[1].degrees())]) + \
                4*sum([get_diff_order(vars_tup[0]), get_diff_order(vars_tup[1])])
    else: 
        return sum(vars_tup[0].degrees()) + 4 * get_diff_order(vars_tup[0])

