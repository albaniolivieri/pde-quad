from sympy import *
from quadratization_copy import is_a_quadratization

# Tests     
u, ux, uxx, uxxx = symbols('u ux uxx uxxx')
w0, w0x, w0xx, w0xxx = symbols('w0 w0x w0xx w0xxx')
u_t, w_0t = symbols('u_t w_0t')

# u_t = u**2 * uxx
# w = u**2
# w_t = 2u * (u**2 * uxx)
V0 = list(map(lambda v, l: (l, poly(v, [u, ux, uxx])), [1, u, ux, uxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx], [1, u, ux, uxx, w0, w0x, w0xx]))
w0t = [(w_0t, poly(2*u**3*uxx, [u, ux, uxx]))]
assert is_a_quadratization(V0, w0t) == [Eq(w_0t, w0*w0xx - w0x**2/2)]

# u_t = u * (ux**2 + u * uxx)
# w = u**2
# w_t = 2u**2 * (ux**2 + u * uxx)
V1 = list(map(lambda v, l: (l, poly(v, [u, ux, uxx])), [1, u, ux, uxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx], [1, u, ux, uxx, w0, w0x, w0xx]))
w0t1 = [(u_t, poly(u*ux**2 + u**2*uxx, [u, ux, uxx])), (w_0t, poly(2*u**2*ux**2 + 2*u**3*uxx, [u, ux, uxx]))]
assert is_a_quadratization(V1, w0t1)#  == [Eq(u_t, u*w0xx/2), Eq(w_0t, w0*w0xx)]

# u_t = u * (3 ux * uxx + u * uxxx + 1)
# w = u**2
# w_t = 2 * u**2 * (2 ux * uxx + u * uxxx + 1)
V2 = list(map(lambda v, l: (l, poly(v, [u, ux, uxx, uxxx])), [1, u, ux, uxx, uxxx, u**2, 2*u*ux, 2*ux**2+2*u*uxx, 6*ux*uxx + 2*u*uxxx],
    [1, u, ux, uxx, uxxx, w0, w0x, w0xx, w0xxx]))
w0t2 = [(u_t, poly(u*(3*ux*uxx + u*uxxx + 1), [u, ux, uxx, uxxx])), (w_0t, poly(2*u**2*(3*ux*uxx + u*uxxx + 1), [u, ux, uxx, uxxx]))]
assert is_a_quadratization(V2, w0t2) #== [Eq(u_t, ux*w0xx/2), Eq(w_0t, w0*w0xxx + 2*w0)]