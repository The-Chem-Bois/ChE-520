import numpy as np
from calculator_functions import get_vapor_pressure


def shortcut_flash(P, T, antoine_coeffs, eps):
    """
    Shortcut method if you do not know the feed compositions! Will return only split fractions.

    P: Pressure in mmHg
    T: Temperature in Kelvin
    antoine_coeffs: Array of antoine coefficients A, B, and C for each key component.
    eps: Split fraction and index of one of the key components.
    """

    p_vaps = get_vapor_pressure(antoine_coeffs, T)

    n = eps.position - 1
    epsilon = eps.value

    k = p_vaps/P

    alpha = k/k[n]

    epsilons = alpha*epsilon/(1+(alpha - 1)*epsilon)

    return (epsilons)