import sympy as sp
import time
import statistics
from qbee import *


if __name__ == '__main__':
    u, ux, uxx, uxxx = functions("u, ux, uxx, uxxx")
    system = [
        (u, ux**3),
        (ux, 3 * ux**2 * uxx),
        (uxx, 6 * ux * uxx**2 + 3 * ux**2 * uxxx)
    ]

    times = []
    for i in range(10):
        ti = time.time()
        quadr_system = polynomialize_and_quadratize(system, input_der_orders={uxxx: 0})
        times.append(time.time()-ti)
        if quadr_system:
            order = quadr_system
            
    print(f'order: {order}', f'avg: {statistics.mean(times)}', f'std dev: {statistics.stdev(times)}')
        
        