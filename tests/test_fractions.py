import sys
sys.path.append("..")

import random
from sympy import apart, symbols, simplify, Add, fraction, expand
from algorithm.fractions import decompose_fraction 

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

dens_symp = [fraction(term)[1] for term in Add.make_args(sympy_decomp)]

expr, rel, q_syms = decompose_fraction(n/d, [w])

print('Fraction decomposition:', expr, '\nRelations:', rel)

final_expr = expr.subs(list(zip(q_syms, [1/fac for fac in rel])))

print('Checking if sympy decomposition is equal to our decomposition:', 
      f'\n{final_expr} - ({sympy_decomp}) =', simplify(expand(final_expr) - expand(sympy_decomp)))

if (simplify(expand(final_expr) - expand(sympy_decomp)) == 0):
    print('Test passed')
else:
    print('Test failed')

