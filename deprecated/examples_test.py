# In this file we test verifying a quadratization for some toy and practical PDE examples 
from sympy import symbols, Function
from sympy import Derivative as D
from deprecated.check_quad_test import test_quad 

# examples

t, x = symbols('t x')
u = Function('u')(t,x)

tests = []

#u_t = u**2 * ux 
#w = u**2
#w_t = 2u * (u**2 * ux + u)
ut1 = u**2*D(u, x) 
w01 = u*D(u, x)
w02 = u**3
tests.append(test_quad([(u, ut1)], [w01, w02], 3))

ut2 = u**2*D(u, x, 2)
w02 = u**2
tests.append(test_quad([(u, ut2)], [w02], 2))

ut2p = D(u, x)**2*u
w02p = u**2
tests.append(test_quad([(u, ut2p)], [w02p], 1))

ut3 = u * (D(u, x)**2 + u * D(u, x, 2))
w03 = u**2
tests.append(test_quad([(u, ut3)], [w03], 2))

ut4 = u * (3 * D(u, x) * D(u, x, 2) + u * D(u, x, 3) + 1)
w04 = u**2
tests.append(test_quad([(u, ut4)], [w04], 3))

#Dym 
ut5 = u**3 * D(u, x, 3)
tests.append(test_quad([(u, ut5)], [u**3, u * D(u, x)**2], 3))

ut6 = D(u, x)**3 + u**3
tests.append(test_quad([(u, ut6)], [u**2, D(u,x)**3], 3))

ut6 = D(u, x)**3 
tests.append(test_quad([(u, ut6)], [D(u, x)**3, D(u,x)*D(u,x,2)], 3))

u1 = Function('u1')(t,x)
u1t = u1**3 * D(u1, x, 1)
ut5 = u**3 * D(u, x, 3)
tests.append(test_quad([(u, ut5), (u1, u1t)], [u**3, u * D(u, x)**2, u1**3], 3))

# Solar wind
v = Function('v')(t,x)
vinv = Function('vinv')(t,x)
vt = 7 * D(v, x) * vinv - 5 * D(v, x)
vinv_t = -vt * vinv**2
tests.append(test_quad([(v, vt), (vinv, vinv_t)], [D(v, x) * vinv, D(v, x) * vinv**2], 3))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))    
