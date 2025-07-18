from __future__ import annotations

from typing import Optional, TYPE_CHECKING
from functools import reduce
import sympy as sp
from sympy import Derivative as D
from sympy.polys.rings import PolyElement

if TYPE_CHECKING:
    from .fraction_decomp import FractionDecomp


def get_order(set_derivs: list[sp.Symbol]) -> int:
    """Gets the maximum order of derivatives in a list of derivatives

    Parameters
    ----------
    set_derivs
        List with the derivatives as expressions

    Returns
    -------
    int
        the maximum order of derivatives
    """

    max_order = 0
    for deriv in set_derivs:
        max_order = reduce(max, [der.args[1][1] for der in deriv.atoms(D)], max_order)
    return max_order


def reduction(pol1: PolyElement, pol2: PolyElement) -> PolyElement:
    """Reduces the first polynomial by the second one

    Parameters
    ----------
    pol1 : PolyElement
        the polynomial to be reduced
    pol2 : PolyElement
        the polynomial to reduce by

    Returns
    -------
    PolyElement
        the reduced polynomial
    """
    coef = pol1[1].coeff_monomial(pol2[2])
    if coef != 0:
        new_pol = pol1[1] - coef * pol2[1]
        return (pol1[0] - coef * pol2[0], new_pol, new_pol.LM())
    else:
        return pol1


def reduction_sparse(
    pol1: tuple[sp.Expr, PolyElement], pol2: tuple[sp.Expr, PolyElement]
) -> tuple[sp.Expr, PolyElement, PolyElement]:
    """Reduces the polynomial pol1 by the polynomial pol2 in the sparse representation

    Parameters
    ----------
    pol1
        the polynomial to be reduced as a tuple with the symbol, the polynomial and the leading monomial
    pol2
        the polynomial to reduce by as a tuple with the symbol, the polynomial and the leading monomial

    Returns
    -------
    tuple[sp.Expr, PolyElement, PolyElement]
        the reduced polynomial as a tuple with the symbol, the polynomial and the leading monomial
    """
    if pol2[1] != 0:
        coef = pol1[1].coeff(pol2[2])
        if coef != 0:
            new_pol = pol1[1] - coef * pol2[1]
            return (
                pol1[0] - coef.as_expr() * pol2[0],
                new_pol,
                new_pol.leading_monom(),
            )
    return pol1


def diff_dict(
    pol: PolyElement,
    dic: dict,
    frac_decomp: Optional[FractionDecomp] = None,
    order: Optional[int] = 1,
) -> PolyElement:
    """differentiates a polynomial by a dictionary of variables

    Parameters
    ----------
    pol
        polynomial to be differentiated
    dic
        dictionary with the variables to differentiate by
    frac_decomp : optional
        the fraction decomposition of the PDE if it has one (default is None)
    order : optional
        the order of differentiation (default is 1)

    Returns
    -------
    PolyElement
        the differentiated polynomial
    """
    deriv = pol
    for _ in range(1, order + 1):
        deriv = sum(deriv.diff(k) * v for (k, v) in dic.items())
    if frac_decomp:
        return pol.ring(frac_decomp.try_reduce(deriv))
    return deriv


def remove_vars(
    list_vars: list[PolyElement], accum_vars: list[PolyElement], axis: int
) -> list[PolyElement]:
    """Removes variables if they are already proposed in the quadratization or if they are of degree less than two

    Parameters
    ----------
    list_vars
        Set of variables to check
    accum_vars
        List with all variables proposed up to this point
    axis
        Axis to check

    Returns
    -------
    list[tuple]
        the filtered list of variables
    """
    for i in range(len(list_vars)):
        if len(list_vars[i]) > 1:
            if list_vars[i][axis] in accum_vars or (
                sum(list_vars[i][axis].degrees()) <= 1
            ):
                list_vars[i] = (list_vars[i][int(not axis)],)
            elif list_vars[i][axis] == list_vars[i][int(not axis)]:
                list_vars[i] = (list_vars[i][axis],)
    return list_vars


def powerset(iterable: list[PolyElement]) -> list[PolyElement]:
    """Returns the powerset of a set

    Parameters
    ----------
    iterable
        Set to get the powerset from

    Returns
    -------
    list[PolyElement]
        the powerset of the set
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)))


def get_diff_order(pol: PolyElement) -> int:
    """Returns the order of the highest derivative in a polynomial

    Parameters
    ----------
    pol
        Polynomial to check

    Returns
    -------
    int
        the order of the highest derivative in the polynomial
    """
    derivs = [x for x in pol.ring.gens if pol.diff(x) != 0]
    order = 0
    for var in derivs:
        if str(var)[-1].isnumeric():
            order += int(str(var)[-1])
    return order
