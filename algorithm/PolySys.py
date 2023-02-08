from .quadratization import is_quadratization
from .utils import diff_dict

class PolySys:
    
    def __init__(self, dic_t, dic_x, pde, order, var_indep, poly_vars, new_vars=[]): 
        self.dic_t = dic_t
        self.dic_x = dic_x
        self.new_vars = new_vars
        self.pde_eq = pde
        self.order = order
        self.var_indep = var_indep
        self.poly_vars = poly_vars
        
    def set_new_vars(self, new_vars):
        self.new_vars = new_vars
        
    def get_quad(self):   
        new_vars_named = [(symbols(f'w_{i}'), pol) for i, pol in enumerate(self.new_vars)] 
        new_vars_t, new_vars_x = self.differentiate_dict(new_vars_named)  
        deriv_t = new_vars_t + self.pde_eq   
        V = [(1, list(self.dic_t.keys())[0].ring(1))] \
            + [(symbols(f'{sym}'), sym) for sym in self.poly_vars] \
            + new_vars_named + new_vars_x
        return is_quadratization(V, deriv_t)
    
    def differentiate_dict(self, named_new_vars):
        deriv_t = []
        deriv_x = []
        
        for name, expr in named_new_vars:
            deriv_t.append((symbols(f'{name}t'), diff_dict(expr, self.dic_t)))
        
        for name, expr in self.new_vars:
            for i in range(1, order + 1):
                deriv_x.append((symbols(f'{name}{self.var_indep}{i}'), 
                                diff_dict(expr, self.dic_x, i)))
                
        return deriv_t, deriv_x