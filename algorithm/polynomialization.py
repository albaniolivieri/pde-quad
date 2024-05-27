from sympy import symbols, Function, exp, sin, cos, Pow, Symbol, Integer
from sympy import Derivative as D

non_pol_funcs = [exp, sin, cos, Pow]

def polynomialize(func_eq, first_indep=symbols('t'), second_indep=symbols('x')):
    # pde_sys = [(lhs, cancel(rhs)) for lhs, rhs in pde_sys]
    count = 0
    j = 0
    new_vars = []
    while j in range(len(func_eq)):
        while True:
            new_var = get_non_pol_var(func_eq[j][1])
            # print('pase la funcion y esta es la new_var', new_var)
            if not new_var:
                break
            new_vars.append((new_var, Function(f'p_{count}')(first_indep, second_indep)))
            for i in range(len(func_eq)):
                func_eq[i] = (func_eq[i][0], func_eq[i][1].subs(new_vars))        
            new_eq = get_new_eq(new_var, func_eq, first_indep).subs(new_vars)
            func_eq.append((Function(f'p_{count}')(first_indep, second_indep), new_eq))
            count += 1
        j += 1
    return func_eq, new_vars
    
def get_non_pol_var(expr, new_var = None):
    args = expr.args
    # print("expr", expr, "args", args, 'new_var', new_var)
    if expr.func in non_pol_funcs:
        if expr.func == Pow:
            if args[1] == int(args[1]):
                return new_var
        for arg in args:
            # print('arg', arg)
            return get_non_pol_var(arg, expr)
    else: 
        if expr.func == Symbol or expr.func == Integer:
            # print('entre aqui en symbol')
            return new_var
        else: 
            for arg in args:
                # print('entre aqui en add o mult')
                new_var = get_non_pol_var(arg, new_var)
                if new_var != None: 
                    return new_var
    return new_var

def get_new_eq(new_var, func_eq, first_indep):
    # print('new_var', new_var)
    refac = [(D(func, first_indep), eq) for func, eq in func_eq]
    wt = D(new_var, first_indep).doit().subs(refac).doit()
    # print('wt', wt)
    return wt

    
    