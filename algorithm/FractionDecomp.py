from sympy import fraction, factor_list, groebner, symbols, sympify
from functools import reduce
from .utils import ring_to_expr, diff_dict
from fractions import Fraction


class FractionDecomp:
    """
    A class to perform a fraction decomposition for a PDE system

    ...

    Attributes
    ----------
    pde : list[tuple]
        A list that represents the PDE system
    rels : list[tuple]
        A list with all the relations introduced by the fraction decomposition
    groeb_rels : list[sympy.Expr]
        A list that represents the ideal I set for the Groebner basis
    q_syms : list[sympy.Symbol]
        A list with all the symbols introduced by the fraction decomposition

    Methods
    -------
    get_frac_decomp(pde_sys, syms)
        Performs the fraction decomposition for the PDE system
    diff_frac(rel, dic, n_diff=1)
        Calculates the derivative of a fraction
    try_reduce(expr)
        Reduces an expression using the Groebner basis        
    """

    def __init__(self, pde_sys, pol_syms, consts):
        """   

        Parameters
        ----------
        pde_sys : list[tuple]
            A list that represents the PDE system
        syms : list[sympy.Symbol]
            A list with all symbols of the PDE system

        Returns
        -------
        None

        """
        new_pde, groeb_rels, rel_list, q_syms = self.get_frac_decomp(
            pde_sys, pol_syms, consts)
        self.pde = new_pde
        self.rels = rel_list
        self.groeb_rels = groeb_rels
        self.q_syms = q_syms

    def get_frac_decomp(self, pde_sys, pol_syms, consts):
        """
        Returns the reduced PDE system, the Groebner basis, the fraction relations and 
        the symbols introduced by the fraction decomposition.

        Parameters
        ----------
        pde_sys : list[tuple]
            A list that represents the PDE system
        syms : list[sympy.Symbol]
            A list with all symbols of the PDE system

        """
        q_symb, rel_list, coef_den = [], {}, []
        i = 0
        for k in range(len(pde_sys)):
            n, d = fraction(pde_sys[k][1])
            coef_den.append(1)
            d_factor = factor_list(d, gens=pol_syms)
            for x in n.as_coefficients_dict().values():
                if x.is_Float:
                    rat_coef = Fraction(str(x)).denominator
                    n *= rat_coef
                    coef_den[k] = coef_den[k]*rat_coef
            coef_den[k] = coef_den[k]*d_factor[0]
            pde_sys[k] = (pde_sys[k][0], n)
            if d != 1:
                for j in range(len(d_factor[1])):
                    rel = d_factor[1][j][0]
                    if not rel_list or rel not in rel_list.values():
                        q = symbols(f'q_{i}')
                        rel_list[q] = rel
                        q_symb.append(q)
                        i += 1
                    else:
                        key_list = list(rel_list.keys())
                        val_list = list(rel_list.values())
                        q = key_list[val_list.index(rel)]
                    pde_sys[k] = (pde_sys[k][0], pde_sys[k][1] * q**d_factor[1][j][1])
        groeb_rels = [rel * rel_list[rel] - 1 for rel in rel_list]
        if groeb_rels:
            groeb_base = groebner(groeb_rels, pol_syms+q_symb+consts, order='lex')
            for k in range(len(pde_sys)):
                pde_sys[k] = (pde_sys[k][0],
                              groeb_base.reduce(pde_sys[k][1])[1] / coef_den[k])
            groeb_rels = [rel.as_expr() for rel in groeb_base._basis]
        rel_list = list(zip(rel_list.keys(), rel_list.values()))
        return pde_sys, groeb_rels, rel_list, q_symb
    
    def rels_as_poly(self, R):
        """
        Transforms the relations to polynomials

        Parameters
        ----------
        R : sympy.PolyRing
            The polynomial ring
        """
        self.groeb_rels = [R.ring_new(rel) for rel in self.groeb_rels]

    def diff_frac(self, rel, dic, n_diff=1):
        """
        Returns the derivative of a fraction

        Parameters
        ----------
        rel : tuple
            A tuple with the relation to be differentiated
        dic : dict
            A dictionary with the differentiation rules
        consts : list[sympy.Symbol]
            A list with all the symbol constants of the PDE system
        """
        q, den = rel
        deriv_var = den.ring(q)
        deriv_num, deriv_den = den.ring(1), den
        for _ in range(1, n_diff + 1):
            deriv_num = (diff_dict(deriv_num, dic) * deriv_den -
                         diff_dict(deriv_den, dic) * deriv_num)
            deriv_den = den**2
            deriv_var = deriv_var**2

        return den.ring(self.try_reduce(deriv_num*deriv_var))

    def try_reduce(self, expr):
        """
        Returns the reduced expression using the Groebner basis

        Parameters
        ----------
        expr : sympy.Expr
            The expression to be reduced
        """
        if not self.groeb_rels:
            return expr
        return expr.div(self.groeb_rels)[1] 
