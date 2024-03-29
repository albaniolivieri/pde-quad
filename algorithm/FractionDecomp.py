from sympy import fraction, factor_list, groebner, symbols, xring, FractionField, QQ
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
    groeb_base : sympy.GroebnerBasis
        The Groebner basis calculated with respect to the relations
    q_syms : list[sympy.Symbol]
        A list with all the symbols introduced by the fraction decomposition

    Methods
    -------
    get_frac_decomp(pde_sys, syms)
        Performs the fraction decomposition for the PDE system
    diff_frac(rel, dic, consts, n_diff=1)
        Calculates the derivative of a fraction
    try_reduce(expr, consts)
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
        new_pde, groeb_base, rel_list, q_syms = self.get_frac_decomp(
            pde_sys, pol_syms, consts)
        self.pde = new_pde
        self.rels = rel_list
        self.groeb_base = groeb_base
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
        q_symb, rel_list, coef_den = [], [], []
        i = 0
        for k in range(len(pde_sys)):
            n, d = fraction(pde_sys[k][1])
            if d == 1:
                continue
            factor_exp, q_exp = [], []
            coef_den_k = 1
            d_factor = factor_list(d, gens=pol_syms)
            for x in n.as_coefficients_dict().values():
                if x.is_Float:
                    rat_coef = Fraction(str(x)).denominator
                    n *= rat_coef
                    coef_den_k *= rat_coef
            coef_den_k *= d_factor[0]
            for j in range(len(d_factor[1])):
                # Gleb: are you sure you really want this story with `coef`?
                # It seems that you could have ignired it, it would just got to the denominators
                # of the coefficients
                factor_exp.append((d_factor[1][j][0], d_factor[1][j][1]))
                rel = d_factor[1][j][0]
                if rel not in rel_list:
                    rel_list.append((symbols(f'q_{i}'), rel))
                    q_symb.append(symbols(f'q_{i}'))
                    i += 1
            # Gleb: It seems that the next double for-loop and reduce aim at writing 
            # the denominator in terms q's. Couldn't you just do this right away during the above
            # loop over the factors? I think this would be simpler
            for q, expr in rel_list:
                for factor, exp in factor_exp:
                    if expr == factor:
                        q_exp.append((q, exp))
            coef_den.append(coef_den_k)
            pde_sys[k] = (pde_sys[k][0], reduce(
                lambda y, x: y * x[0]**x[1], q_exp, 1) * n)
        groeb_rels = [q * fac - 1 for q, fac in rel_list]
        if groeb_rels:
            groeb_base = groebner(groeb_rels, pol_syms+q_symb+consts, order='lex')
            for k in range(len(pde_sys)):
                pde_sys[k] = (pde_sys[k][0],
                              groeb_base.reduce(pde_sys[k][1])[1] / coef_den[k])
        else:
            groeb_base = None
        return pde_sys, groeb_base, rel_list, q_symb

    def diff_frac(self, rel, dic, consts, n_diff=1):
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
        q = den.ring(q)
        deriv_num, deriv_den = den.ring(1), den
        for _ in range(1, n_diff + 1):
            deriv_num = (diff_dict(deriv_num, dic) * deriv_den -
                         diff_dict(deriv_den, dic) * deriv_num)
            deriv_den = den**2
            deriv_var = q**2
        return den.ring(self.try_reduce(deriv_num*deriv_var, consts))

    def try_reduce(self, expr, consts):
        """
        Returns the reduced expression using the Groebner basis

        Parameters
        ----------
        expr : sympy.Expr
            The expression to be reduced
        consts : list[sympy.Symbol]
            A list with all the symbol constants of the PDE system
        """
        # Gleb: perhaps we have discussed this, but do we really need to have a conversion to expression here?
        coefs_denom = 1
        for x in expr.coeffs():
            coefs_denom *= fraction(x.as_expr())[1]
        return self.groeb_base.reduce(
            ring_to_expr(expr*coefs_denom, self.groeb_base.gens, consts))[1]/coefs_denom
