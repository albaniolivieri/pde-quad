from sympy import Eq, pprint
from .utils import reduction_sparse


def is_quadratization(V, deriv, frac_decomp):
    """Verifies if all variables in V are a quadratization for the system in deriv.
    Also prints the quadratization

    Parameters
    ----------
    V : list[sympy.PolyElement]
        List V with all variables/symbols
    deriv : list[sympy.PolyElement]
        List with all equations of the PDE system
    frac_decomp : FractionDecomp
        Fraction decomposition of the PDE system

    Returns
    -------
    tuple
        a tuple with a boolean value and the quadratization found (if not found, returns the NS set)
    """
    if frac_decomp.groeb_rels:
        V2 = list(set((m1[0] * m2[0], m1[1].ring(frac_decomp.try_reduce(
            m1[1] * m2[1]))) for m1 in V for m2 in V))
    else:
        V2 = list(set((m1[0] * m2[0], m1[1] * m2[1]) for m1 in V for m2 in V))

    V2_poly, names = [], []
    for name, polyn in V2:
        names.append(name)
        V2_poly.append(polyn)
    quad, NS = [], []

    V2_red = reduce_set(V2)

    for name, pol in deriv:
        if pol not in V2_poly:
            result = is_linear_combination(V2_red, pol, name)
            if type(result) == tuple:
                NS.append((name, result[1][1]))
            else:
                quad.append(Eq(name, result))
        else:
            quad.append(Eq(name, names[V2_poly.index(pol)]))
    if NS != []:
        # for printing problematic monomials, uncomment next line:
        # for i in range(len(NS)): pprint(f'NS for expr {NS[i][0]}: {NS[i][1]}')
        return (False, NS)

    print("\nQuadratization:")
    for exprs in quad:
        pprint(exprs)
    return (True, quad)


def reduce_set(V2):
    """Reduces the V^2 set following the Gauss elimination method 

    Parameters
    ----------
    V2 : list[sympy.PolyElement]
        List with V^2 set

    Returns
    -------
    list
        a list with the reduced V^2 set 
    """
    for i in range(len(V2)):
        for j in range(i):
            V2[i] = reduction_sparse(V2[i], V2[j])
        if V2[i][1] != 0:
            lead_coeff = V2[i][1].coeff(V2[i][1].leading_monom())
            V2[i] = (V2[i][0] / lead_coeff.as_expr(), 
                     V2[i][1] * (1 / lead_coeff), 
                     V2[i][1].leading_monom())
            for j in range(i):
                V2[j] = reduction_sparse(V2[j], V2[i])
    return [(a[0], a[1]) for a in V2]


def is_linear_combination(V2, der_pol, name):
    """Checks if a certain polynomial is a linear combination of the set V^2

    Parameters
    ----------
    V2 : list[sympy.PolyElement]
        List with V^2 set
    der_pol : sympy.PolyElement
        Polynomial to check

    Returns
    -------
    tuple or sympy.PolyElement
        if der_pol is a linear combination of V^2, returns the resulting polynomial from done operations.
        if it is not, returns a tuple with the boolean False and a tuple that represents der_pol polynomial
    """
    der_tuple = (0, der_pol, der_pol.leading_monom())
    V2 = [(name, pol, pol.leading_monom()) for name, pol in V2]
    for i in range(len(V2)):
        der_tuple = reduction_sparse(der_tuple, V2[i])
        if der_tuple[1] == 0:
            return -der_tuple[0]
    return (False, der_tuple)