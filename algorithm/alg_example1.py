from sympy import symbols, Poly, prod, Eq, solve

u, ux, uxx, w_0xx = symbols('u ux uxx w_0xx')

p = Poly(2*ux**2+2*u*uxx)

monom = [prod(x**k for x, k in zip(p.gens, mon)) for mon in p.monoms()]

coeffs = p.coeffs()

#print(monom)
#print(coeffs)

V = [1, u, uxx, ux, u**2, p, 2*u*ux] # u*ux, ux**2, u*uxx]
V2 = []

for var in V: 
    for var2 in V:
        mult = var*var2
        #print(mult)
        if mult not in V2:
            V2.append(mult)
            print(f"{mult} \n")

#print("V2:", V2.join\n)
der_t = u**3*uxx

NS = []

if der_t not in V2: 
    NS.append(der_t)

print(NS)

#dic_poly = {"first" : solve(Eq(w_0xx, p), 2*ux**2), "second" : solve(Eq(w_0xx, p), 2*u*uxx)}

#print(dic_poly)
