import sys
from sympy import symbols, simplify, expand
from sympy import Derivative as D

sys.path.append("..")

from algorithm.utils import get_order
from algorithm.check_manual_quad import test_quadratization

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
        quad_vars.extend([(symbols(f'q_{i}{var_indep}{j}'), D(frac_vars[i], var_indep, j).doit())
                          for j in range(1, n + 1)])
    return quad_vars


def test_quad(func_eq: list, new_vars: list, n_diff: int, frac_vars: list = [], vars_from_alg: bool = False):
    """Test the proposed quadratization of a given PDE

    Parameters
    ----------
    func_eq : list
        Tuples with the symbol and equations of the PDE
    new_vars : list
        List of proposed new variables
    n_diff : int
        The number of second variable differentiations to do
    frac_vars : list, optional
        List of fraction decomposition variables

    Returns
    -------
    bool
        True if the quadratization is correct, False otherwise
    """
    x_var = [symbol for symbol in func_eq[0]
             [0].free_symbols if symbol != symbols('t')].pop()
    
    max_order = get_order(list(zip(*func_eq))[1])
    
    refac = []
    for fun, _ in func_eq:
        refac += [(symbols(f'{fun.name}_{x_var}{i}'), D(fun, x_var, i))
                  for i in range(n_diff + max_order + 1, 0, -1)] + [(symbols(fun.name), fun)]
        
    new_vars_expr = new_vars
    if vars_from_alg:
        quad_prop = [expr.subs(refac) for expr in new_vars]
        frac_vars = [(q, 1/expr.subs(refac)) for q, expr in frac_vars]
        new_vars = [expr.subs(frac_vars) for expr in quad_prop]
    else: 
        frac_refac = [(expr, q) for q, expr in frac_vars]
        input_refac = [(expr, name) for name, expr in refac]
        new_vars_expr = [expr.subs(frac_refac + input_refac) for expr in new_vars]
        
    var_dic = [(symbols(f'w_{i}'), new_vars[i]) for i in range(len(new_vars))]
    total_vars = (new_vars, [rel for _, rel in frac_vars])
    quad_vars = differentiate_x(x_var, total_vars, n_diff)
    deriv_t = differentiate_t(func_eq, [(var, expr.subs(frac_vars)) for var, expr in var_dic] + frac_vars) \
        + [(symbols(eqs[0].name + '_t'), eqs[1]) for eqs in func_eq]
    refac += quad_vars + frac_vars
    exprs_orig = [expr for _, expr in deriv_t]
    
    results = test_quadratization(func_eq, new_vars_expr, n_diff)

    if not results[0] and not results[1]:
        print("\nQuadratization not found")
        return False
    
    for i in range(len(exprs_orig)):
        print('Checking equation:', results[1][i])
        if simplify(exprs_orig[i] - results[1][i].rhs.subs(refac)) != 0: 
            print('Test failed: expressions are not equal')
            print('Equation: ', results[1][i])
            print('Original expression: ', expand(simplify(exprs_orig[i])))
            print('Quad expression: ', expand(simplify(results[1][i].rhs.subs(refac))))
            print('Substraction: ', simplify(exprs_orig[i] - results[1][i].rhs.subs(refac)))
            return False

    print('Test passed')
    return True
