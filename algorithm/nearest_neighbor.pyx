from queue import PriorityQueue
from collections import deque
from .var_selection import prop_new_vars
from .utils import shrink_quad

def nearest_neighbor(poly_syst, sort_fun, new_vars=[]):
    pq = PriorityQueue()
    
    pq.put((len(new_vars), 0, new_vars))
    NS_queue = deque()
    node_count = 0
    count = 0
    
    while not pq.empty():
        # print('size', pq.qsize()) 
        new_vars = pq.get()[2]
        poly_syst.set_new_vars(new_vars)
        result_quad = poly_syst.try_make_quadratic()
        node_count += 1
    
        if result_quad[0]:
            shrinked_quad = shrink_quad(new_vars, poly_syst)
            return shrinked_quad, node_count
        else:
            NS_queue.append((new_vars, result_quad[1]))
            if pq.qsize() <= 1: 
                new_vars, NS = NS_queue.popleft()
                prop_vars = prop_new_vars(NS, new_vars, sort_fun)
                for p_vars in prop_vars: 
                    pq.put((len(new_vars + list(p_vars)), count, new_vars + list(p_vars)))
                    count += 1