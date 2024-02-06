from sympy import xring, symbols, QQ
from sympy import Derivative as D
from .quadratization import is_quadratization
from .utils import diff_dict, get_order, diff_frac

class PolySys: 
    """
    A class used to represent a PDE system as polynomial expressions

    ...

    Attributes
    ----------
    dic_t : dict
        a dictionary that stores time derivatives derivatives
    dic_x : dict
        a dictionary that stores spatial derivatives
    new_vars : list[sympy.Symbol]
        a list of new variables introduced for the PDE quadratization
    pde_eq : list[tuple]
        a list of tuples containing the PDE equations in the form
        (<left hand side as symbol>, <right hand side as polynomial>)
    order : int
        an integer representing the order of second variable derivatives
    var_indep : sympy.Symbol
        a Symbol that represents the second independent variable
    poly_vars : sympy.PolyElement
        all symbol names as polynomial ring
    expr_frac : list[tuple]
        a list of tuples containing the fraction variables and their expressions
    frac_ders : list[sympy.PolyElement]
        a list of the t derivatives of the fraction variables

    Methods
    -------
    build_ring(pde_sys, n_diff, var_indep, new_vars=[])
        Builds the polynomial ring for all the expressions
    get_dics(func_eq, symb, eqs_pol, order, max_order)
        Builds dictionary that links a symbol with 
        its symbol derivative
    set_new_vars(new_vars)
        Sets the attribute new_vars to parameter passed 
    try_make_quadratic()
        Gets the quadratization of the PDE
    differentiate_dict(named_new_vars)
        Builds two dictionaries with new variables derivatives in first
        and second variable
    """ 
    
    def __init__(self, pde_sys, n_diff, var_indep, new_vars=[], vars_frac=[]): 
        """
        Parameters
        ----------
        pde_sys : list[tuple]
            Tuples with the symbol and expression of PDE
        n_diff : int
            The number of second variable differentiations to do 
        var_indep : sympy.Symbol
            The symbol of the second independent variable
        new_vars : list, optional
            List of proposed new variables
        vars_frac : list, optional
            List of new variables of the fraction decomposition
        """
        max_order = get_order([expr for _, expr in pde_sys])
        new_vars_pol = new_vars
       
        poly_syms, eqs_pol, new_vars_pol, expr_frac = self.build_ring(pde_sys, n_diff, var_indep, max_order, vars_frac, new_vars)
            
        dic_t, dic_x, frac_ders = self.get_dics(pde_sys, poly_syms, eqs_pol, n_diff, max_order, expr_frac)
        
        self.dic_t = dic_t
        self.dic_x = dic_x
        self.new_vars = new_vars_pol
        self.pde_eq = eqs_pol
        self.order = n_diff
        self.var_indep = var_indep
        self.poly_vars = poly_syms
        self.expr_frac = expr_frac
        self.frac_ders = frac_ders
    
    
    def build_ring(self, func_eq, order, var_indep, max_order, frac_vars, new_vars=None):
        """Returns equation symbols and expressions expressed as polynomials. 
        If new_vars parameter is passed, it also returns the new variables as polynomials.

        Parameters
        ----------
        func_eq : list[tuple]
            Tuples with the symbol and expression of PDE
        order : int
            The number of second variable differentiations to do 
        var_indep : sympy.Symbol
            The symbol of the second independent variable
        max_order : int
            Max order of derivatives in all the system
        new_vars : list, optional
            List of proposed new variables
        """
        
        # der_subs is used for derivatives and functions substitution to sympy symbols
        der_subs = []
        for fun, _ in func_eq:
            der_subs += [(D(fun, var_indep, i), symbols(f'{fun.name}_{var_indep}{i}')) 
                    for i in range(max_order, 0, -1)] + [(fun, symbols(fun.name))] 
        
        poly_vars = []
        for fun, _ in func_eq:
            poly_vars.append(f'{fun.name}')
            poly_vars.extend([f'{fun.name}_{var_indep}{i}' for i in range(1, max_order + order + 1)])
            
        for var, _ in frac_vars:
            poly_vars.append(f'{var.name}') 
        
        R, pol_sym = xring(poly_vars, QQ)
        
        frac_subs = [(1/rel, var) for var, rel in frac_vars]
        
        func_eq_frac = [(fun, eq.subs(frac_subs)) for fun, eq in func_eq]   

        expr_pol = [(symbols(f'{fun.name}_t'), R.ring_new(eq.subs(der_subs))) for fun, eq in func_eq_frac] 
        
        expr_frac = [(symbols(f'{var.name}'), R.ring_new(rel.subs(der_subs))) for var, rel in frac_vars]    
        
        # if the new variables are passed as sympy expressions  
        # they are also added to the polynomial ring
        vars_pol = []
        if new_vars != None: 
            vars_pol = [R.ring_new(new_vars[i].subs(der_subs)) for i in range(len(new_vars))]
        
        return pol_sym, expr_pol, vars_pol, expr_frac
    
    def get_dics(self, func_eq, symb, eqs_pol, order, max_order, expr_frac): 
        """Returns a tuple with two dictionaries that map the equations symbols 
        with their respective derivatives and a list with the derivatives of the
        fraction variables

        Parameters
        ----------
        func_eq : list[tuple]
            Tuples with the symbol and expression of PDE
        symb : sympy.PolyElement
            PDE symbols in polynomial ring 
        eqs_pol : list[tuple]
            Tuples containing the PDE equations in the form
        new_vars : list, optional
            List of proposed new variables 
        order : int
            The number of second variable differentiations to do 
        """
        dic_x = {}
        dic_t = {}
        
        der_order = max_order + order
        for j in range(len(func_eq)):
            for i in range(der_order): 
                dic_x[symb[i + (der_order+1)*j]] = symb[i + (der_order+1)*j + 1]
                last = i + (der_order+1)*j
        
        count = last
        for i in range(len(expr_frac)):
            dic_x[symb[count+2]] = diff_frac(symb[0].ring(1), expr_frac[i][1], symb[last+2], dic_x)
            count += 1
        
        for k in range(len(func_eq)):
            for i in range(der_order - 1): 
                if i != 0:
                    dic_t[symb[i + (der_order+1)*k]] = diff_dict(dic_t[symb[i - 1 + (der_order+1)*k]], dic_x)
                else: 
                    dic_t[symb[(der_order+1)*k]] = eqs_pol[k][1]
        
        frac_ders = []            
        for i in range(len(expr_frac)):
            frac_der_t = diff_frac(symb[0].ring(1), expr_frac[i][1], symb[last+2], dic_t)
            frac_ders.append(frac_der_t)
            dic_t[symb[last+2]] = frac_der_t
            last += 1
            
        return dic_t, dic_x, frac_ders
        
    def set_new_vars(self, new_vars):
        """Returns none as it only sets with a new value the new_vars attribute

        Parameters
        ----------
        new_vars : list
            List of proposed new variables 
        """
        self.new_vars = new_vars
        
    def try_make_quadratic(self):  
        """Returns a tuple with a bool and the quadratization if it is a quadratization. If not,
        returns the expressions that could not be quadratized.
        """
        new_vars_named = [(symbols(f'w_{i}'), pol) for i, pol in enumerate(self.new_vars)] 
        new_vars_t, new_vars_x = self.differentiate_dict(new_vars_named) 
        frac_der_t = [(symbols(f'q_{i}t'), pol) for i, pol in enumerate(self.frac_ders)]
        deriv_t = new_vars_t + self.pde_eq + frac_der_t
        V = [(1, list(self.dic_t.keys())[0].ring(1))] \
            + [(symbols(f'{sym}'), sym) for sym in self.poly_vars] \
            + new_vars_named + new_vars_x
        return is_quadratization(V, deriv_t, self.expr_frac, self.poly_vars)
   
    def differentiate_dict(self, named_new_vars):
        """Returns two dictionaries that map the new variables with their respective
        derivatives in the first and second variable 

        Parameters
        ----------
        named_new_vars : list[tuple]
            List of proposed new variables with their respective symbol
        """
        deriv_t = []
        deriv_x = []
        
        for name, expr in named_new_vars:
            deriv_t.append((symbols(f'{name}t'), diff_dict(expr, self.dic_t)))
        
        for name, expr in named_new_vars:
            for i in range(1, self.order + 1):
                deriv_x.append((symbols(f'{name}{self.var_indep}{i}'), 
                                diff_dict(expr, self.dic_x, i)))
                
        return deriv_t, deriv_x
