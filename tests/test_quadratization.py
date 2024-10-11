import unittest
import sys
import math
from sympy import symbols, simplify, expand, Function 
from sympy import Derivative as D

sys.path.append("..")

from algorithm.quadratize import quadratize
from algorithm.utils import ring_to_expr, get_order
from algorithm.check_manual_quad import check_quadratization
from algorithm.var_selection import by_order_degree, by_fun, by_degree_order

class TestCase():
    
    def __init__(self, func_eq, n_diff, max_der_order=None) -> None:
        self.func_eq = func_eq
        self.n_diff = n_diff
        self.max_der_order = max_der_order
        
class TestQuadratization(unittest.TestCase):
    
    def setUp(self):
        self.t, self.x = symbols('t x')
        self.u = Function('u')(self.t, self.x)
        self.v = Function('v')(self.t, self.x)
        self.omega = symbols('omega', constant=True)
        
        self.test_cases_quad = []
        # Example: u_t = u**2 * u_xx + 2, v_t = D(v, x, 2)/u**3 + u
        self.test_cases_quad.append(TestCase([(self.u, - D(self.v, self.x, 1) * self.u**3 - 2 * D(self.v, self.x, 1) * self.u**2), 
                  (self.v, D(self.v, self.x, 1) * self.u - 2 * D(self.v, self.x, 1))], 2, max_der_order=2))
        self.test_cases_quad.append(TestCase([(self.u, self.u**3 * D(self.u, self.x, 3))], 3))
        self.test_cases_quad.append(TestCase([(self.u, D(self.u, self.x)**3 + self.u**3)], 2))
        self.test_cases_quad.append(TestCase([(self.u, D(self.u, self.x)**4)], 3, max_der_order=2))
        self.test_cases_quad.append(TestCase([(self.u, D(self.u, self.x)**3 * self.u)], 2, max_der_order=3))
        self.test_cases_quad.append(TestCase([(self.u, D(self.u, self.x)**3)], 2, max_der_order=3))
        
        self.test_cases_rat = []
        self.test_cases_rat.append(TestCase([(self.u, 3.4 * D(self.u, self.x) * D(self.v, self.x)), 
                  (self.v, D(self.v, self.x) / self.v - round((math.pi), 5) * D(self.v, self.x))], 4))
        self.test_cases_rat.append(TestCase([(self.u, self.u**2 * D(self.u, self.x, 2) + 2), 
                  (self.v, D(self.v, self.x, 2)/self.u**3 + self.u)], 3))
        # self.test_cases_rat.append(TestCase([(self.u, 1 / (5 * (self.u + 1)))], 1))
        # self.test_cases_rat.append(TestCase([(self.u, 1 / (0.6 * self.u + 1.3) ** 2)], 3))
        self.test_cases_rat.append(TestCase([(self.u, 1 / (self.u**2 + 1))], 3))
        self.test_cases_rat.append(TestCase([(self.u, D(self.u, self.x) / (self.u + 1))], 3))
        self.test_cases_rat.append(TestCase([(self.u, 1 / (self.u + 1) ** 2 + 1 / (self.u - 1))], 4, max_der_order=4))
        self.test_cases_rat.append(TestCase([(self.u, 1 / (self.u**2) - 0.5 * self.u + 1)], 4, max_der_order=4))
        self.test_cases_rat.append(TestCase([(self.u,  -D(self.u, self.x, 2) / self.u - self.u**2 - self.u + 5), 
                  (self.v, self.u / self.v - self.v + 5)], 4))
        
        self.test_cases_coeff_sym = []
        self.test_cases_coeff_sym.append(TestCase([(self.u, self.omega * self.u**3 * D(self.u, self.x, 3))], 3))
        self.test_cases_coeff_sym.append(TestCase([(self.u, 1 / (self.omega * (self.u + 1)))], 1))
         
    def transform_new_vars(self, new_vars, frac_vars):
        quad_prop = list(map(lambda x: ring_to_expr(x), new_vars))
        frac_vars = list(map(lambda x: ring_to_expr(x), frac_vars))
        return quad_prop, frac_vars
    
    def differentiate_t(self, funcs_eqs, new_vars):
        deriv_t = []
        refac = [(D(deriv[0], symbols('t')), deriv[1]) for deriv in funcs_eqs]
        for i in range(len(new_vars)):
            wt = D(new_vars[i][1], symbols('t')).doit().subs(refac)
            deriv_t.append((symbols(f'{new_vars[i][0]}_t'), wt.doit()))
        return deriv_t
    
    def differentiate_x(self, new_vars, n_diff):
        vars_prop, frac_vars = new_vars
        quad_vars = []
        for i in range(len(vars_prop)):
            quad_vars.extend([(symbols(f'w_{i}{self.x}{j}'), D(vars_prop[i], self.x, j).doit())
                            for j in range(1, n_diff + 1)] + [(symbols(f'w_{i}'), vars_prop[i])])
        for i in range(len(frac_vars)):
            quad_vars.extend([(symbols(f'q_{i}{self.x}{j}'), D(frac_vars[i], self.x, j).doit())
                            for j in range(1, n_diff + 1)])
        return quad_vars
    
    def refac_expr(self, test_case, new_vars, frac_vars):
        max_order = get_order(list(zip(*test_case.func_eq))[1])
        refac = []
        for fun, _ in test_case.func_eq:
            refac += [(symbols(f'{fun.name}_{self.x}{i}'), D(fun, self.x, i))
                    for i in range(test_case.n_diff + max_order + 1, 0, -1)] + [(symbols(fun.name), fun)]
            
        quad_prop = [expr.subs(refac) for expr in new_vars]
        frac_vars = [(q, 1/expr.subs(refac)) for q, expr in frac_vars]
        new_vars = [expr.subs(frac_vars) for expr in quad_prop]
        
        return new_vars, frac_vars, refac
    
    def construct_quadratic_PDE(self, new_vars, frac_vars, test_case, refac):
        var_dic = [(symbols(f'w_{i}'), new_vars[i]) for i in range(len(new_vars))]
        total_vars = (new_vars, [rel for _, rel in frac_vars])
        quad_vars = self.differentiate_x(total_vars, test_case.n_diff)
        deriv_t = self.differentiate_t(test_case.func_eq, [(var, expr.subs(frac_vars)) for var, expr in var_dic] + frac_vars) \
            + [(symbols(eqs[0].name + '_t'), eqs[1]) for eqs in test_case.func_eq]
        refac += quad_vars + frac_vars
        exprs_orig = [expr for _, expr in deriv_t]
        return refac, exprs_orig
         
    def quadratization_test(self, search_alg, test_cases, sort_heur=by_fun):
        for test in test_cases:
            print('\nTest case: ', test.func_eq)
            quad = quadratize(test.func_eq, test.n_diff, sort_fun=sort_heur, search_alg=search_alg, max_der_order=test.max_der_order)
            self.assertIsNotNone(quad, 'Quadratization not found')
            quad_prop, frac_vars, _ = quad
            self.assertTrue(quad_prop or frac_vars, 'Quadratization not found')
            
            quad_prop_expr, frac_vars_expr = self.transform_new_vars(quad_prop, frac_vars)
            new_vars, frac_vars, refac = self.refac_expr(test, quad_prop_expr, frac_vars_expr)
            refac_new_vars, exprs_orig = self.construct_quadratic_PDE(new_vars, frac_vars, test, refac)
            
            results = check_quadratization(test.func_eq, quad_prop, test.n_diff)
            for i in range(len(exprs_orig)):
                rewritten_result = results[1][i].rhs.subs(refac_new_vars)
                self.assertEqual(simplify(exprs_orig[i] - rewritten_result), 0, 
                                f"Test failed: expressions are not equal \n" + \
                                f"Equation: {results[1][i]} \n" + \
                                f"Original expression: {expand(simplify(exprs_orig[i]))} \n" + \
                                f"Quad expression: {expand(simplify(rewritten_result))} \n" + \
                                f"Substraction: {simplify(exprs_orig[i] - rewritten_result)}") 
                
    def test_branch_and_bound(self):
        self.quadratization_test(search_alg='bnb', test_cases=self.test_cases_quad)
   
    def test_nearest_neighbor(self):
        self.quadratization_test(search_alg='nn', test_cases=self.test_cases_quad)
    
    def test_rational_pdes(self):
        self.quadratization_test(search_alg='nn', test_cases=self.test_cases_rat)
    
    def test_symbolic_coeff(self):
        self.quadratization_test(search_alg='nn', test_cases=self.test_cases_coeff_sym)
        
    
if __name__ == '__main__':
    unittest.main()