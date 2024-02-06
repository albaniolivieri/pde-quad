#!/usr/bin/python
import time
import sys

from sympy import *
from sympy import Derivative as D

sys.path.append("..")
from algorithm import check_manual_quad as sparse
from deprecated import check_quad_gauss as gauss
from deprecated import check_quad_matrix as matrix

t, x = symbols('t x')
u = Function('u')(t,x)
u1 = Function('u1')(t,x)

u1t = u1**3 * D(u1, x, 1)
ut5 = u**3 * D(u, x, 3)

if len(sys.argv) > 1:
    if sys.argv[1] == 'sparse':
        for i in range(1, 11):
            ti = time.time()
            sparse.test_quadratization([(u, ut5), (u1, u1t)], [u**3, u * D(u, x)**2, u1**3], i)
            print("Diff order:", i, "\nTime:", time.time() - ti)
    elif sys.argv[1] == 'gauss':
        for i in range(1, 11):
            ti = time.time()
            gauss.get_quadratization([(u, ut5), (u1, u1t)], [u**3, u * D(u, x)**2, u1**3], i)
            print("Diff order:", i, "\nTime:", time.time() - ti)
    elif sys.argv[1] == 'matrix':
        for i in range(1, 11):
            ti = time.time()
            matrix.get_quadratization([(u, ut5), (u1, u1t)], [u**3, u * D(u, x)**2, u1**3], i)
            print("Diff order:", i, "\nTime:", time.time() - ti)
    
