import sys
sys.path.append("..")

import random
from sympy import apart, symbols, simplify, Add, fraction, expand
from algorithm.FractionDecomp import FractionDecomp

w = symbols('w')

def generate_rand_poly(sym, degree):
    p = 0
    for i in range(degree+1):
        p += random.randint(0,degree) * sym**i
    return p 

n = generate_rand_poly(w, 3)
d = generate_rand_poly(w, 3)

print('Numerator generated:', n)
print('Denominator generated:', d)

sympy_decomp = apart(n/d)

frac_decomp = FractionDecomp([(symbols('u'), n/d)], [symbols('u'), w])

print('Fraction decomposition:', frac_decomp.pde[0][1], '\nRelations:', frac_decomp.rels)

final_expr = frac_decomp.pde[0][1].subs([(q, 1/fac) for q, fac in frac_decomp.rels])

print('final_expr', final_expr)
print('sympy_decomp', sympy_decomp)

print('Checking if sympy decomposition is equal to our decomposition:', 
      f'\n{final_expr} - ({sympy_decomp}) =', simplify(expand(final_expr) - expand(sympy_decomp)))

if (simplify(expand(final_expr) - expand(sympy_decomp)) == 0):
    print('Test passed')
else:
    print('Test failed')

