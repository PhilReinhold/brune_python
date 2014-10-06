import operator
from numpy import linspace, array, argmin
from numpy.polynomial.polynomial import Polynomial, polydiv, polymul, polysub
from scipy.optimize import fmin
from sympy import symbols, poly
from vectfit import vectfit_auto

def threshold_arr(arr, threshold=1e-4):
    arr = array(arr)
    abs_arr = abs(arr)
    abs_thresh = threshold * max(abs_arr)
    ret = arr.copy()
    ret[abs_arr < abs_thresh] = 0.0
    return ret

def get_rational_rep(y, s, n_poles, threshold=1e-4):
    poles, residues, d, h = vectfit_auto(y, s, n_poles=n_poles, show=True)
    x = symbols('x')
    denominator = reduce(operator.mul, [x-p for j, p in enumerate(poles)])
    numerator = 0
    for i, r in enumerate(residues):
        term = r
        for j, p in enumerate(poles):
            if i != j:
                term *= x - p
        numerator += term

    numerator += d*denominator + h*x*denominator

    num_coefs = threshold_arr(map(complex,poly(numerator, x).all_coeffs()), threshold)
    denom_coefs = threshold_arr(map(complex,poly(denominator, x).all_coeffs()), threshold)

    numerator = Polynomial(list(reversed(num_coefs)))
    denominator = Polynomial(list(reversed(denom_coefs)))
    print numerator
    print '-'*30
    print denominator
    print numerator.degree(), denominator.degree()
    return numerator, denominator

def get_local_minima(num, denom):
    z = lambda x: num(x)/denom(x)
    d_num = num.deriv()
    d_denom = denom.deriv()
    dzdx_num = polysub(polymul(d_num, denom) - polymul(num, d_denom))
    roots = dzdx_num.roots()
    local_minima = [z(x) for x in roots]
    min_i = argmin(local_minima)
    min_z = local_minima[min_i]
    if z(0) < min_z:
        return 0
    roots[argmin(local_minima)]

def test_get_rational_rep():
    from matplotlib import pyplot as plt
    f = linspace(0, 10, 501)[1:]
    s = 1j*f
    y = sum([1/(r + 1j*zc*(f/f0 - f0/f)) for r, zc, f0 in [(.2, 10, 3), (.5, 17, 6)]])
    plt.figure(0)
    numerator, denominator = get_rational_rep(y, s, 2)
    n1 = polymul(numerator.coef, [0, 1])
    print polydiv(denominator.coef, n1)
    y2 = numerator(s) / denominator(s)

    plt.figure(1)
    plt.plot(f, y.imag)
    plt.plot(f, y2.imag)
    plt.figure(2)
    plt.plot(f, y.real)
    plt.plot(f, y2.real)
    plt.show()

def single_step_polydiv(a, b):
    'a should be of larger degree than b'
    b2 = polymul(b, [0, 1])
    r, new_a2 = polydiv(a, b2)
    new_a = polymul(new_a2, [0, 1])
    return r, new_a

def reduce_order(num, denom, s):
    if num.degree() > denom.degree():
        # pole at infinity
        r, new_num = single_step_polydiv(num.coef, denom.coef)
        return ('Pole', 'Inf', r), Polynomial(new_num), denom

    if denom.degree() > num.degree():
        # zero at infinity
        r, new_denom = single_step_polydiv(denom.coef, num.coef)
        return ('Zero', 'Inf', r), num, Polynomial(new_denom)

    if denom.coef[0] == 0:
        # pole at zero
        num_inverse = list(reversed(num.coef))
        denom_inverse = list(reversed(num.coef))
        r, new_num_inverse = single_step_polydiv(num_inverse, denom_inverse)
        new_num = list(reversed(new_num_inverse))
        return ('Pole', 'Zero', r), Polynomial(new_num), denom

    if num.coef[0] == 0:
        # zero at zero
        num_inverse = list(reversed(num.coef))
        denom_inverse = list(reversed(num.coef))
        r, new_denom_inverse = single_step_polydiv(denom_inverse, num_inverse)
        new_denom = list(reversed(new_denom_inverse))
        return ('Zero', 'Zero', r), num, Polynomial(new_denom)


    z = lambda x: num(x)/denom(x)
    s1 = get_min(num, denom)
    w1 = s1.imag
    z1 = z(s1)
    r = z1.real
    l1 = z1.imag / w1
    # Remove the resulting zero at w1
    np = polysub(num.coef, polymul(denom.coef, [r, 1j*l1]))
    err_k = lambda k: polydiv(denom - polymul([0, 2*k], np))[0]
    fmin()


    dnum = num.deriv()
    ddenom = denom.deriv()
    dzds = lambda x: (dnum(x)*denom(x) - num(x)*ddenom(x))/(denom(x)**2)
    min_idx = argmin(z(s))
    s1 = s[min_idx]
    w1 = s1.imag
    z1 = z(s1)
    r = z1.real
    l1 = z1.imag / w1
    zc = w1 * dzds(s1).imag / 2.
    l2 = zc / w1
    c2 = 1 / (zc * w1)
    l3 = -l1*l2 / (l1+l2)
    return ('Brune', r, l1, l2, c2, l3),



if __name__ == '__main__':
    test_get_rational_rep()
