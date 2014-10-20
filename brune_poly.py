"""
brune.py -- Implementation of Brune network synthesis algorithm
"""
import matplotlib
matplotlib.use('TKAgg')
import operator
import numpy as np
from numpy.polynomial.polynomial import Polynomial

def product(xs):
    return reduce(operator.mul, xs, 1)

def poles_to_rational_rep(poles, residues, d, h):
    """
    Take the parameters produced by vectfit and construct a representation
    of the form f(s) / g(s), returning f and g
    """
    pole_ps = [Polynomial([-p, 1]) for p in poles]
    denom = product(pole_ps)
    num = sum(r * denom / p for r, p in zip(residues, pole_ps))
    num += Polynomial([d, h]) * denom
    # Since we're a PR function, should be safe to convert these to real values
    num = Polynomial(num.coef.real)
    denom = Polynomial(denom.coef.real)
    return num, denom

def get_real_polynomial(num, denom):
    """get a rational representation of Re[z(iw)]"""
    # 1.A.i: Write num(s) as num_w(w) == num(iw), similarly for denom
    N = num.degree() + 1
    num_w_coefs = np.array([1j**n for n in range(N)])*num.coef
    denom_w_coefs = np.array([1j**n for n in range(N)])*denom.coef

    # 1.A.ii: write Re[z[(iw)]] as
    # (Re[num]Re[denom] + Im[num]Im[denom])/(Re[denom]**2 + Im[denom]**2)
    re_num_w, im_num_w = map(Polynomial, (num_w_coefs.real, num_w_coefs.imag))
    re_denom_w, im_denom_w = map(Polynomial, (denom_w_coefs.real, denom_w_coefs.imag))
    re_z_w_num = re_num_w*re_denom_w + im_num_w*im_denom_w
    re_z_w_denom = re_denom_w**2 + im_denom_w**2
    return re_z_w_num, re_z_w_denom

def get_polynomial_minimum(num, denom):
    d_num = num.deriv()
    d_denom = denom.deriv()
    w_mins = (d_num*denom - num*d_denom).roots().real
    p_w_mins = w_mins[w_mins > 0]
    i0 = np.argmin(num(p_w_mins) / denom(p_w_mins))
    return p_w_mins[i0]

def brune_extraction(num, denom):
    assert num.degree() == denom.degree()
    # Identify s at which re[z(iw)] is minimized
    re_z_w_num, re_z_w_denom = get_real_polynomial(num, denom)
    w0 = get_polynomial_minimum(re_z_w_num, re_z_w_denom)
    s0 = 1j*w0
    z0 = num(s0) / denom(s0)

    # Extract resistor and inductor, producing lossless pole
    r = z0.real
    l1 = z0.imag / w0
    num -= denom * Polynomial([r, l1])

    # Remove the resulting zero
    zero_p = Polynomial([w0**2, 0, 1])
    new_num = num / zero_p
    s_new_num = Polynomial([0, 1])*new_num
    s_new_num_r = s_new_num % zero_p
    denom_r = denom % zero_p
    res = denom_r / s_new_num_r
    new_denom = (denom - res*s_new_num) / zero_p

    yc = res.coef[0] / w0
    c2 = yc / w0
    l2 = 1 / (yc * w0)
    l3 = -l1 * l2 / (l1 + l2)

    final_num = new_num - new_denom*Polynomial([0, l3])

    assert final_num.degree() == new_denom.degree()
    return (r, l1, c2, l2, l3), final_num, new_denom

def brune_stage(s, params, next_z):
    r, l1, c2, l2, l3 = params
    def par(x, y):
        return 1/(1/x + 1/y)
    return r + s*l1 + par(s*l2 + 1/(s*c2), s*l3 + next_z)

if __name__ == '__main__':
    import vectfit
    import matplotlib.pyplot as plt

    poles = [-1 + 100j, -1 - 100j, -2 + 500j, -2 - 500j]
    residues = [14 + 3j, 14 - 3j, 7 + 1j, 7 - 1j]
    offset = 0.5
    w = np.linspace(1, 700, 1000)
    s = 1j*w
    z = vectfit.model(s, poles, residues, offset, 0)
    num, denom = poles_to_rational_rep(poles, residues, offset, 0)
    print num, denom
    plt.plot(w, z.real)
    plt.plot(w, (num(s)/denom(s)).real)
    r_num, r_denom = get_real_polynomial(num, denom)
    plt.plot(w, r_num(w)/r_denom(w))
    plt.yscale('log')
    params, new_num, new_denom = brune_extraction(num, denom)
    brune_z = brune_stage(s, params, new_num(s)/new_denom(s))
    plt.plot(w, brune_z.real)
    plt.show()
