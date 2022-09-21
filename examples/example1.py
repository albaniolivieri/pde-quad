import sympy as sp
from qbee import *


if __name__ == '__main__':
    u, ux, uxx, uxxx = functions("u, ux, uxx, uxxx")
    system = [
        (u, ux * u**2),
        (ux, uxx * u**2 + 2 * ux**2 * u),
        (uxx, uxxx * u**2 + 6 * u * ux * uxx + 2 * ux**3)
    ]

    quadr_system = polynomialize_and_quadratize(system, input_der_orders={uxxx: 0})
    if quadr_system:
        print(quadr_system)