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
    flag = False
    quad_temp = None
    
    while not pq.empty():
        # print('size', pq.qsize()) 
        new_vars = pq.get()[2]
        if quad_temp:
            if len(new_vars) >= len(quad_temp):
                continue
        poly_syst.set_new_vars(new_vars)
        result_quad = poly_syst.try_make_quadratic()
        node_count += 1
    
        if result_quad[0]:
            shrinked_quad = shrink_quad(new_vars, poly_syst)
            if flag: 
                if quad_temp:
                    if len(shrinked_quad) < len(quad_temp):
                        quad_temp = shrinked_quad
            else:
                quad_temp = shrinked_quad
            flag = True
            while len(NS_queue) > 0:
                new_vars_ns, NS = NS_queue.popleft()
                if len(new_vars_ns) + 1 < len(quad_temp):
                    prop_vars = prop_new_vars(NS, new_vars_ns, sort_fun)
                    for p_vars in prop_vars: 
                        if len(new_vars_ns + list(p_vars)) < len(quad_temp):
                            pq.put((len(new_vars_ns + list(p_vars)), count, new_vars_ns + list(p_vars)))
                            count += 1
        else:
            if not flag:
                NS_queue.append((new_vars, result_quad[1]))
                if pq.qsize() <= 1:  
                        new_vars, NS = NS_queue.popleft()
                        prop_vars = prop_new_vars(NS, new_vars, sort_fun)
                        for p_vars in prop_vars: 
                            pq.put((len(new_vars + list(p_vars)), count, new_vars + list(p_vars)))
                            count += 1
    return quad_temp, node_count
