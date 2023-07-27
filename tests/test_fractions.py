import sys
sys.path.append("..")

import random
from sympy import ring, QQ, apart, symbols, simplify, Add, fraction, apart_list
from algorithm.fractions import decompose_fraction 
from algorithm.utils import ring_to_expr, expr_to_ring

w = symbols('w')

def generate_rand_poly(sym, degree):
    p = 0
    # change it to generate factors of some degree for not having irrepla
    for i in range(degree+1):
        p += random.randint(0,degree) * sym**i
    return p 

n = generate_rand_poly(w, 3)
d = generate_rand_poly(w, 3)

print('numerator', n)
print('denominator', d)

sympy_decomp = apart(n/d)
list_apart = apart_list(n/d)
print('sym_decomp', sympy_decomp)
print('lis_apart', list_apart)

dens_symp = [fraction(term)[1] for term in Add.make_args(sympy_decomp)]

expr, rel, q_syms = decompose_fraction(n/d, [w])

print('myDdecomp', expr, rel)

final_expr = expr.subs(list(zip(q_syms, [1/fac for fac in rel])))

#print(True) if rel in dens_symp else print(False)
print(f'{final_expr} - ({sympy_decomp}) =', simplify(final_expr - sympy_decomp))

#together: for combining denominators and then do the comparison for every factor  
# i dont check exactly that decomposition is the same
# then check the lists