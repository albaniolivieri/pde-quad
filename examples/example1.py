import sympy as sp
import statistics
import time
from qbee import *


if __name__ == '__main__':
    u, ux, uxx, uxxx = functions("u, ux, uxx, uxxx")
    system = [
        (u, ux * u**2),
        (ux, uxx * u**2 + 2 * ux**2 * u),
    ]

    times = []
    for i in range(10):
        ti = time.time()
        quadr_system = polynomialize_and_quadratize(system, input_der_orders={uxx: 0})
        times.append(time.time()-ti)
        if quadr_system:
            order = quadr_system
            
    print(f'order: {order}', f'avg: {statistics.mean(times)}', f'std dev: {statistics.stdev(times)}')