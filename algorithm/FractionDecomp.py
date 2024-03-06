from sympy import fraction, factor_list, groebner, symbols
from functools import reduce
from .utils import ring_to_expr, diff_dict

class FractionDecomp:
    def __init__(self, pde_sys, syms):        
        new_pde, groeb_base, rel_list, q_syms = self.get_frac_decomp(pde_sys, syms)
        self.pde = new_pde
        self.rels = rel_list
        self.groeb_base = groeb_base
        self.q_syms = q_syms
        
    def get_frac_decomp(self, pde_sys, syms):
        q_symb, rel_list = [], []
        i = 0
        for k in range(len(pde_sys)):
            n, d = fraction(pde_sys[k][1])
            if d == 1:
                continue
            factor_exp = []
            q_exp = []
            d_factor = factor_list(d)
            coef = d_factor[0]
            for j in range(len(d_factor[1])):
                if coef != 1:
                    factor_exp.append((d_factor[1][j][0]*coef, d_factor[1][j][1]))
                    rel = d_factor[1][j][0]*coef
                    coef = 1
                else:
                    factor_exp.append((d_factor[1][j][0], d_factor[1][j][1]))
                    rel = d_factor[1][j][0]
                if rel not in rel_list:
                    rel_list.append((symbols(f'q_{i}'), rel))
                    q_symb.append(symbols(f'q_{i}'))
                    i += 1
            for q, expr in rel_list:
                for factor, exp in factor_exp:
                    if expr == factor:
                        q_exp.append((q, exp))
            pde_sys[k] = (pde_sys[k][0], reduce(lambda y, x: y * x[0]**x[1], q_exp, 1) * n)
        groeb_rels = [q * fac - 1 for q, fac in rel_list]
        if groeb_rels:
            groeb_base = groebner(groeb_rels, q_symb + syms, order='lex')
            new_pde = list(map(lambda x: (x[0], groeb_base.reduce(x[1])[1]), pde_sys))
        else: 
            groeb_base = None
            new_pde = pde_sys
        return new_pde, groeb_base, rel_list, q_symb
        
    def diff_frac(self, rel, dic, consts, n_diff = 1):
        q, den = rel
        q = den.ring(q)
        deriv_num, deriv_den = den.ring(1), den
        for _ in range(1, n_diff + 1):
            deriv_num = (diff_dict(deriv_num, dic) * deriv_den - diff_dict(deriv_den, dic) * deriv_num) 
            deriv_den = den**2
            deriv_var = q**2
        return den.ring(self.try_reduce(deriv_num*deriv_var, consts))
        
    def try_reduce(self, expr, consts):
        return self.groeb_base.reduce(ring_to_expr(expr, self.groeb_base.gens, consts))[1]
        
                
        
        