import sympy as sp
from qbee import *


if __name__ == '__main__':
    u, ux, uxx, uxxx = functions("u, ux, uxx, uxxx")
    system = [
        (u, ux * u**2 + u),
        (ux, 2 * u * ux**2 + u**2 * uxx + ux),
        (uxx, 2 * ux**3 + 6 * ux * uxx * u + u**2 * uxxx + uxx)
    ]

    quadr_system = polynomialize_and_quadratize(system, input_der_orders={uxxx: 0})
    if quadr_system:
        print(quadr_system)