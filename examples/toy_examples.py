from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize

#tests
t, x = symbols('t x')
u = Function('u')(t,x)

orders = []

# u_t = u**2 * ux 
# w = u**2
# w_t = 2u * (u**2 * ux + u)
ut1 = u**2 * D(u, x)

# u_t = u**2 * ux + u**3
ut2 = u**2 * D(u, x) + u**3

# u_t = ux**2 * u 
ut3 = u * D(u, x)**2

# u_t = ux**3 
ut4 = D(u, x)**3

# u_t = u**2 * uxx
# w = u**2
# w_t = 2u * (u**2 * uxx)
ut5 = u**2 * D(u, x, 2)

# u_t = u * (ux**2 + u * uxx)
# w = u**2
# w_t = 2u**2 * (ux**2 + u * uxx)
ut6 = u * (D(u, x)**2 + u * D(u, x, 2))

# u_t = u * (2 ux * uxx + u * uxxx + 1)
# w = u**2
# w_t = 2 * u**2 * (2 ux * uxx + u * uxxx + 1)
ut7 = u * (3 * D(u, x) * D(u, x, 2) + u * D(u, x, 3) + 1)

pdes = [ut1, ut2, ut3, ut4, ut5, ut6, ut7]

for eq in pdes: 
    times = []
    for i in range(10):
        ti = time.time()
        order = quadratize([(u, eq)], 3)[1]
        times.append(time.time() - ti)
        
    print(f'result for {eq}', f'order: {order}', f'avg: {statistics.mean(times)}',
        f'std dev: {statistics.stdev(times)}')