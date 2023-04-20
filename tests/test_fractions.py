import sys
sys.path.append("..")

import random
from sympy import ring, QQ, apart, symbols, simplify
from algorithm.fractions import decompose_fraction 
from algorithm.utils import ring_to_expr, expr_to_ring

# R, x, y, z = ring('x y z', QQ)
# n = (x**2 + y *z - z**5)
# d = ((x + y) * (y + z)**2 * (z + x))
# n_expr, expr_syms = ring_to_expr([x,y,z], n)
# d_expr, expr_syms = ring_to_expr([x,y,z], d)
# print(decompose_fraction(n_expr/d_expr, expr_syms))

w = symbols('w')

def generate_rand_poly(sym, degree):
    p = 0
    for i in range(degree+1):
        p += random.randint(0,degree) * sym**random.randint(0,degree)
    return p 

n = generate_rand_poly(w, 3)
d = generate_rand_poly(w, 3)

print('numerator', n)
print('denominator', d)

sympy_decomp = apart(n/d)

expr, rel, q_syms = decompose_fraction(n/d, [w])

final_expr = expr.subs(list(zip(q_syms, [1/fac for fac in rel])))

print(f'{final_expr} - ({sympy_decomp}) =', simplify(final_expr - sympy_decomp))