from sympy import xring, symbols, QQ, expand
from sympy import Derivative as D
from .quadratization import is_quadratization
from .utils import diff_dict, get_order, diff_frac
from .fractions import get_frac_decomp

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
    
    def __init__(self, pde_sys, n_diff, vars_indep, new_vars=[], frac_vars=[]): 
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
            List of new variables from the fraction decomposition
        """
        max_order = get_order([expr for _, expr in pde_sys])
        
        self.first_indep, self.sec_indep = vars_indep
        self.constants = []
        self.frac_vars = frac_vars
       
        poly_syms, eqs_pol, new_vars_pol, frac_decomps = self.build_ring(pde_sys, n_diff, vars_indep, max_order, new_vars)
            
        dic_t, dic_x, frac_der_t = self.get_dics(pde_sys, poly_syms, eqs_pol, n_diff, max_order, frac_decomps)
        
        self.dic_t = dic_t
        self.dic_x = dic_x
        self.new_vars = new_vars_pol
        self.pde_eq = eqs_pol
        self.order = n_diff
        self.poly_vars = poly_syms
        self.expr_frac = frac_decomps
        self.frac_der_t = frac_der_t
    
    def build_ring(self, func_eq, order, var_indep, max_order, new_vars=None):
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
            der_subs += [(D(fun, var_indep[1], i), symbols(f'{fun.name}_{var_indep[1].name}{i}')) 
                    for i in range(max_order, 0, -1)] + [(fun, symbols(fun.name))] 

        symbols_list = [symbol for _, expr in func_eq for symbol in expr.free_symbols]
        
        poly_vars = []
        constants = [symbols(x.name, constant = True) for x in 
                     set(filter(lambda x: (x != self.first_indep) and (x != self.sec_indep), 
                                 symbols_list))]
        poly_vars.extend(constants)
        self.constants = constants
        
        for fun, _ in func_eq:
            poly_vars.append(symbols(fun.name))
            poly_vars.extend([symbols(f'{fun.name}_{var_indep[1].name}{i}') 
                              for i in range(1, max_order + order + 1)])
        
        func_eq = [(lhs, rhs.subs(der_subs)) for lhs, rhs in func_eq] 
        
        frac_decomps = get_frac_decomp(func_eq, poly_vars)
        
        if frac_decomps:
            for decomp in frac_decomps:
                for i in range(len(func_eq)):
                    if func_eq[i][0] in decomp[0]:
                        func_eq[i] = decomp[0]
        
        for decomp in frac_decomps:
            for j in range(len(decomp[2])):
                poly_vars.append(symbols(decomp[2][j].name))
                poly_vars.extend([symbols(f'q_{j}{var_indep[1].name}{i}') 
                                  for i in range(1, max_order + order + 1)])
        
        R, pol_sym = xring(poly_vars, QQ) 

        expr_pol = [(symbols(f'{fun.name}_{var_indep[0]}'), R.ring_new(eq)) for fun, eq in func_eq]
        
        vars_pol = {'frac_vars': [], 'new_vars': []}
        for i in range(len(frac_decomps)):
            frac_decomps[i][1] = [R.ring_new(expr) for expr in decomp[1]]
            frac_decomps[i][2] = [R.ring_new(var) for var in decomp[2]]
            vars_pol['frac_vars'].extend(frac_decomps[i][1])   
        
        # if the new variables are passed as sympy expressions  
        # they are also added to the polynomial ring
        if new_vars != None: 
            new_vars = [var.subs(der_subs) for var in new_vars]
            for decomp in frac_decomps:
                subs_fracs = [(self.frac_vars[i].subs(der_subs), symbols(f'q_{i}')) 
                              for i in range(len(self.frac_vars))]
                #print('subs_fracs', subs_fracs)
                new_vars = [var.subs(subs_fracs) for var in new_vars]
                #print('new_vars', new_vars)
                #new_vars = [decomp[3].reduce(var)[1] for var in new_vars]
            vars_pol['new_vars'] = [R.ring_new(var) for var in new_vars]
        
        return pol_sym, expr_pol, vars_pol, frac_decomps
    
    def get_dics(self, func_eq, symb, eqs_pol, order, max_order, frac_decomps): 
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
        
        constants = [symb[0].ring(const) for const in self.constants]
        symb = list(filter(lambda x: x not in constants, symb))
        
        der_order = max_order + order
        for i in range(len(func_eq)):
            for j in range(der_order): 
                dic_x[symb[j + (der_order+1)*i]] = symb[j + (der_order+1)*i + 1]
                last = j + (der_order+1)*i
        
        frac_ders = []
        count = last+2
        for decomp in frac_decomps:
            j=count
            for i in range(len(decomp[2])):
                dic_x[symb[j + (der_order+1)*i]] = diff_frac(symb[0].ring(1), decomp[1][i], dic_x, decomp[3], decomp[2][i])
                count = count + (der_order+1)*i
                
        for k in range(len(func_eq)):
            for i in range(der_order - 1): 
                if i != 0:
                    dic_t[symb[i + (der_order+1)*k]] = diff_dict(dic_t[symb[i - 1 + (der_order+1)*k]], dic_x)
                else: 
                    dic_t[symb[(der_order+1)*k]] = eqs_pol[k][1]
                    
        count = last+2
        for decomp in frac_decomps:
            for i in range(len(decomp[2])):
                frac_der_t = diff_frac(symb[0].ring(1), decomp[1][i], dic_t, decomp[3], decomp[2][i], constants=self.constants)
                frac_ders.append((symbols(f'{decomp[2][i]}{self.first_indep}'), frac_der_t))
                dic_t[symb[count]] = frac_der_t
                count += der_order+1

        return dic_t, dic_x, frac_ders
        
    def set_new_vars(self, new_vars):
        """Returns none as it only sets with a new value the new_vars attribute

        Parameters
        ----------
        new_vars : list
            List of proposed new variables 
        """
        self.new_vars['new_vars'] = new_vars
        
    def get_frac_vars(self):
        """Returns none as it only sets with a new value the new_vars attribute

        Parameters
        ----------
        new_vars : list
            List of proposed new variables 
        """
        return self.new_vars['frac_vars']
   
    def differentiate_dict(self, named_new_vars):
        """Returns two lists that map the new variables with their respective
        derivatives in the first and second independent variable 

        Parameters
        ----------
        named_new_vars : list[tuple]
            List of proposed new variables with their respective symbol
        """
        deriv_t = []
        deriv_x = []
        
        for name, expr in named_new_vars:
            deriv_t.append((symbols(f'{name}{self.first_indep}'), diff_dict(expr, self.dic_t)))
        
        for name, expr in named_new_vars:
            for i in range(1, self.order + 1):
                deriv_x.append((symbols(f'{name}{self.sec_indep}{i}'), 
                                diff_dict(expr, self.dic_x, i)))
                
        for decomp in self.expr_frac:
            for i in range(len(decomp[2])):
                for j in range(1, self.order + 1):
                    deriv_x.append((symbols(f'{decomp[2][i]}{self.sec_indep}{j}'), 
                                    diff_frac(decomp[1][i].ring(1), decomp[1][i], self.dic_x, 
                                              decomp[3], decomp[2][i], j, constants=self.constants)))               
        return deriv_t, deriv_x
    
    def try_make_quadratic(self):  
        """Returns a tuple with a bool and the quadratization if it is a quadratization. If not,
        returns the expressions that could not be quadratized.
        """
        new_vars_named = [(symbols(f'w_{i}'), pol) for i, pol in enumerate(self.new_vars['new_vars'])] 
        new_vars_t, new_vars_x = self.differentiate_dict(new_vars_named) 
        deriv_t = new_vars_t + self.frac_der_t + self.pde_eq 
        V = [(1, list(self.dic_t.keys())[0].ring(1))] \
            + [(symbols(f'{sym}'), sym) for sym in self.poly_vars] \
            + new_vars_named + new_vars_x 
            
        return is_quadratization(V, deriv_t)
