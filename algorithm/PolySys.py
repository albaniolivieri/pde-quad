from sympy import xring, symbols, QQ, cancel, FractionField, simplify
from sympy import Derivative as D
from .quadratization import is_quadratization
from .utils import diff_dict, get_order
from .FractionDecomp import FractionDecomp

# Gleb: perhaps, we should think about naming at some point:
# this is no longer really a `PolySys` but a rational system instead
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
    # Gleb: shouldn't these just go to dic_t ?
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

    def __init__(self, pde_sys, n_diff, vars_indep, new_vars=None):
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
        # Gleb: What this `vars_frac` are used for ?
        vars_frac : list, optional
            List of new variables from the fraction decomposition
        """
        self.max_order = get_order([expr for _, expr in pde_sys])

        self.first_indep, self.sec_indep = vars_indep
        self.consts = []
        self.order = n_diff

        poly_syms, eqs_pol, new_vars_pol, frac_decomps = self.build_ring(
            pde_sys, new_vars)

        self.frac_decomps = frac_decomps
        
        self.poly_vars = poly_syms
        self.pde_eq = eqs_pol
        self.new_vars = new_vars_pol

        dic_t, dic_x, frac_der_t = self.get_dics(pde_sys)

        self.dic_t = dic_t
        self.dic_x = dic_x
        self.frac_der_t = frac_der_t

    def build_ring(self, func_eq, new_vars):
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
            der_subs += [(D(fun, self.sec_indep, i), symbols(f'{fun.name}_{self.sec_indep.name}{i}'))
                         for i in range(self.max_order, 0, -1)] + [(fun, symbols(fun.name))]

        symbols_list = [symbol for _,
                        expr in func_eq for symbol in expr.free_symbols]
        constants = [symbols(x.name, constant=True) for x in
                     set(filter(lambda x: (x != self.first_indep) and (x != self.sec_indep),
                                symbols_list))]
        self.consts = constants

        poly_vars = []

        for fun, _ in func_eq:
            poly_vars.append(symbols(fun.name))
            poly_vars.extend([symbols(f'{fun.name}_{self.sec_indep.name}{i}')
                              for i in range(1, self.max_order + self.order + 1)])

        # here we substite derivatives symbs and if there is an expr with a rational function,
        # we convert it to the form p/q
        func_eq = [(lhs, cancel(rhs.subs(der_subs))) for lhs, rhs in func_eq]

        frac_decomp = FractionDecomp(func_eq, poly_vars, constants)

        if frac_decomp:
            func_eq = frac_decomp.pde
            for i in range(len(frac_decomp.q_syms)):
                poly_vars.append(frac_decomp.q_syms[i])
                poly_vars.extend([symbols(f'q_{i}{self.sec_indep.name}{k}')
                                  for k in range(1, self.max_order + self.order + 1)])
                
        # print('frac_decomp', frac_decomp.rels)

        QQc = FractionField(QQ, constants)
        R, pol_sym = xring(poly_vars, QQc)
        
        frac_decomp.rels_as_poly(R)

        pde_pol = [(symbols(f'{fun.name}_{self.first_indep}'),
                     R.ring_new(simplify(eq))) for fun, eq in func_eq]
        
        vars_pol = {'frac_vars': [], 'new_vars': []}

        for i in range(len(frac_decomp.rels)):
            frac_decomp.rels[i] = (
                frac_decomp.rels[i][0], R.ring_new(frac_decomp.rels[i][1]))
            vars_pol['frac_vars'].append(frac_decomp.rels[i])

        # if the new variables are passed as sympy expressions
        # they are also added to the polynomial ring
        if new_vars:
            vars_pol['new_vars'] = [R.ring_new(var) for var in new_vars]

        return pol_sym, pde_pol, vars_pol, frac_decomp

    def get_dics(self, func_eq):
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

        der_order = self.max_order + self.order 
        for i in range(len(func_eq)):
            for j in range(der_order):
                dic_x[self.poly_vars[j + (der_order + 1) * i]
                      ] = self.poly_vars[j + (der_order + 1) * i + 1]
                last = j + (der_order + 1) * i

        frac_ders = []
        count = last+2
        rels = self.frac_decomps.rels
        for i in range(len(rels)):
            dic_x[self.poly_vars[count]] = self.frac_decomps.diff_frac(
                rels[i], dic_x) 
            count += (der_order + 1)

        for k in range(len(func_eq)):
            for i in range(der_order):
                if i != 0:
                    dic_t[self.poly_vars[i + (der_order+1)*k]] = diff_dict(
                        dic_t[self.poly_vars[i - 1 + (der_order+1)*k]], dic_x, self.frac_decomps)
                else:
                    dic_t[self.poly_vars[(der_order + 1)*k]
                          ] = self.pde_eq[k][1]

        count = last+2
        for rel in self.frac_decomps.rels:
            frac_der_t = self.frac_decomps.diff_frac(rel, dic_t)
            frac_ders.append(
                (symbols(rel[0].name+self.first_indep.name), frac_der_t))
            dic_t[self.poly_vars[count]] = frac_der_t
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
        """Returns the fraction variables introduced in the system

        Returns
        -------
        list
            a list of the fraction variables introduced
        """
        return self.new_vars['frac_vars']

    def get_max_order(self):
        """Returns the max derivative order of the system

        Returns
        -------
        int
            the max derivative order of the system 
        """
        return self.max_order

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
            deriv_t.append(
                (symbols(f'{name}{self.first_indep}'), diff_dict(expr, self.dic_t, self.frac_decomps)))
            
        for name, expr in named_new_vars:
            for i in range(1, self.order + 1):
                deriv_x.append((symbols(f'{name}{self.sec_indep}{i}'),
                                diff_dict(expr, self.dic_x, order=i, frac_decomp=self.frac_decomps)))

        for rel in self.frac_decomps.rels:
            for j in range(1, self.order + 1):
                deriv_x.append((symbols(f'{rel[0].name}{self.sec_indep}{j}'),
                                self.frac_decomps.diff_frac(
                                    rel, self.dic_x, n_diff=j)))

        return deriv_t, deriv_x

    def try_make_quadratic(self):
        """Returns a tuple with a bool and the quadratization if it is a quadratization. If not,
        returns the expressions that could not be quadratized.
        """
        new_vars_named = [(symbols(f'w_{i}'), pol)
                          for i, pol in enumerate(self.new_vars['new_vars'])]
        new_vars_t, new_vars_x = self.differentiate_dict(new_vars_named)
        deriv_t = new_vars_t + self.frac_der_t + self.pde_eq
        poly_vars = list(filter(lambda x: str(x)[0] != 'q', self.poly_vars))
        
        # print('new_vars', new_vars_named)

        V = [(1, self.poly_vars[0].ring(1))] + [(symbols(f'{sym}'), sym) for sym in poly_vars] \
            + [(q, self.poly_vars[0].ring(q)) for q, _ in self.frac_decomps.rels] \
            + new_vars_named + new_vars_x

        return is_quadratization(V, deriv_t, self.frac_decomps)
