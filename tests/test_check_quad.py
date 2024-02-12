from sympy import symbols, Function, simplify, expand
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm import check_manual_quad as quad
from algorithm.utils import get_order

def differentiate_t(funcs_eqs, new_vars):
    """Differentiate the functions with respect to t
    
    Parameters
    ----------
    funcs_eqs : list[tuple]
        Tuples with the symbol and expression of PDE
    new_vars : list
        List of proposed new variables
        
    Returns
    -------
    list[tuple]
        the differentiated functions with respect to t
    """
    deriv_t = []
    refac = [(D(deriv[0], symbols('t')), deriv[1]) for deriv in funcs_eqs]
    for i in range(len(new_vars)):
        wt = D(new_vars[i][1], symbols('t')).doit().subs(refac)
        deriv_t.append((symbols(f'{new_vars[i][0]}_t'), wt.doit()))
    return deriv_t

def differentiate_x(var_indep, new_vars, n):
    """Differentiates functions with respect to the x variable.
    
    Parameters
    var_indep : sympy.Symbol
        The symbol of the second independent variable
    new_vars : tuple
        Tuple with a list of proposed new variables and 
        another list with fraction relations
    n : int
        The number of second variable differentiations to do
    
    Returns
    -------
    list[tuple]
        the differentiated functions with respect to x
    """
    vars_prop, frac_vars = new_vars
    quad_vars = []
    for i in range(len(vars_prop)):
        quad_vars.extend([(symbols(f'w_{i}{var_indep}{j}'), D(vars_prop[i], var_indep, j).doit()) 
                          for j in range(1, n + 1)] + [(symbols(f'w_{i}'), vars_prop[i])])
    for i in range(len(frac_vars)):
        quad_vars.extend([(symbols(f'q_{i}{var_indep}{1}'), D(frac_vars[i], var_indep, 1).doit())]
                          + [(symbols(f'q{i}'), frac_vars[i])])
    return quad_vars

def test_quad(func_eq, new_vars: list, n_diff: int, frac_vars: list = []):
    """Test the proposed quadratization of a given PDE
    
    Parameters
    ----------
    func_eq : list
        Tuples with the symbol and equations of the PDE
    new_vars : list
        List of proposed new variables
    n_diff : int
        The number of second variable differentiations to do
        
    Returns
    -------
    bool
        True if the quadratization is correct, False otherwise
    """
    x_var = [symbol for symbol in func_eq[0][0].free_symbols if symbol != symbols('t')].pop()
    new_vars = [expr.subs([(symbols(f'q{i}'), frac_vars[i]) for i in range(len(frac_vars))]) 
                for expr in new_vars]
    var_dic = [(symbols(f'w_{i}'), new_vars[i]) for i in range(len(new_vars))] 
    total_vars = (new_vars, frac_vars)
    quad_vars = differentiate_x(x_var, total_vars, n_diff)
    undef_fun = [symbol for symbol, _ in func_eq]
    deriv_t = differentiate_t(func_eq, var_dic) + [(symbols(eqs[0].name + '_t'), eqs[1]) for eqs in func_eq] 
    max_order = max(get_order([der for _, der in deriv_t]), get_order([der for _, der in quad_vars]))
    
    refac = []
    for fun in undef_fun:
        refac += [(symbols(f'{fun.name}_{x_var}{i}'), D(fun, x_var, i)) 
                  for i in range(max_order, 0, -1)] + [(symbols(fun.name), fun)]
    refac += quad_vars
    exprs_orig = [expr for _, expr in deriv_t]
    results = quad.test_quadratization(func_eq, new_vars, n_diff)
    if not results[0]: 
        print("\nQuadratization not found")
        return False 
    
    for i in range(len(exprs_orig)):
        print('passed eq:', results[1][i])
        if expand(exprs_orig[i]) - expand(results[1][i].rhs.subs(refac)) != 0:
            print('Test failed: expressions are not equal')
            print('Equation: ', results[1][i])
            print('Original expression: ', expand(exprs_orig[i]))
            print('Quad expression: ', expand(results[1][i].rhs.subs(refac)))
            return False
        
    print('Test passed')
    return True

