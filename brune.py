from itertools import combinations
from pylab import *
import warnings, sys
from vectfit import vectfit_auto

def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))
warnings.showwarning = customwarn

seterr(all='ignore')


def brune_step(s, z):
    w = imag(s)
    # Extract resistor (real part)
    i1 = nanargmin(real(z))
    r1 = z[i1].real
    s1 = s[i1]
    w1 = w[i1]
    print "pole at", w1, "resistor", r1
    z1 = z - r1

    if i1 == 0:
        # parallel inductor
        real_slope = ((z[1] - z[0]) / (w[1] - w[0])).real
        r_at_zero = max(z[0].real - real_slope * w[0], 0)
        z1 = z - r_at_zero
        #l2 = (z1[0] / s[0]).real
        l2 = ((z[1] - z[0]) / (s[1] - s[0])).real
        stage = (r_at_zero, 0, l2, float('inf'), 0)
        y1 = 1 / z1
        y3 = y1 - 1 / (l2 * s)
        z4 = 1 / y3
        rf = nanmin(real(z4))
        return stage, z4, rf

    if i1 == len(z) - 1:
        # series capacitor
        y = 1 / z
        dy = y[-1] - y[-2]
        ds = s[-1] - s[-2]
        c2 = (dy / ds).real
        stage = (0, 0, 0, c2, 0)
        y3 = y - c2*s
        z4 = 1 / y3
        rf = nanmin(real(z4))
        quadrant_plot(w, 1/(c2*s), "CapacitorZ")
        return stage, z4, rf

    # Extract inductor (imag part)
    l1 = (z1[i1] / s1).real
    z2 = z1 - l1*s

    # Find pole impedance
    dz = (z2[i1+1] - z2[i1-1]).imag
    dw = w[i1+1] - w[i1-1]
    zc = w1 * dz / (2 * dw)
    print "zc", zc, "dz", dz, "dw", dw
    z_res =  1j * zc * (w/w1 - w1/w)
    # z3 *is* undefined at w1, it just is
    z3 = 1 / (1/z2 - 1/z_res)

    # Extract final (negative) inductor
    l2 = zc / w1
    c2 = 1 / (zc * w1)
    l3 = -l1 * l2 / (l1 + l2)
    z4 = z3 - l3*s

    rf = nanmin(real(z4))

    stage = r1, l1, l2, c2, l3
    # quadrant_plot(w, z2, "z2")
    # quadrant_plot(w, z3, "z3")
    quadrant_plot(w, z_res, "z_res")
    return stage, z4, rf

def par(*args):
    div_mask = True
    for x in args:
        div_mask = np.bitwise_and(x != 0, div_mask)

    for x in args:
        if isinstance(x, np.ndarray):
            res = np.zeros_like(x)
            break
    else:
        raise ValueError

    for x in args:
        if isinstance(x, np.ndarray):
            res[div_mask] += 1/x[div_mask]
        else:
            res[div_mask] += 1/x

    res[div_mask] = 1 / res[div_mask]
    return res

def brune_estimate(s, stages, rf=0):
    z = ones(s.shape) * rf
    for stage in reversed(stages):
        try:
            r1, l1, l2, c2, l3, z_stages, y_stages, mask = stage
            post_stages = True
        except:
            r1, l1, l2, c2, l3 = stage
            post_stages = False

        if c2 == float('inf'):
            z = r1 + s*l1 + par(l2*s, l3*s + z)
        else:
            z = r1 + s*l1 + par(l2*s + 1/(c2*s), l3*s + z)
        if post_stages:
            z[mask] = 1 / ((1/z)[mask] + y_stages)
            z[mask] += z_stages
    return z

def quadrant_plot(w, z, name):
    subplot(2,2,1)
    plot(w, real(z), label="Re[%s]" % name)
    legend()
    subplot(2,2,2)
    plot(w, imag(z), label="Im[%s]" % name)
    legend()
    subplot(2,2,3)
    plot(w, real(1/z), label="Re[1/%s]" % name)
    legend()
    subplot(2,2,4)
    plot(w, imag(1/z), label="Im[1/%s]" % name)
    legend()


def foster_preamble(w, s, n_poles, threshold=1e-4, show=False):
    poles, residues, d, h = vectfit_auto(w, s, n_poles, n_iter=100, show=show)
    thresh = threshold * median(poles.real)
    stages = 0
    new_w = w.copy()
    count = 0
    for pole, res in zip(poles, residues):
        if pole.real > thresh:
            print 'lossless pole', res, pole
            count += 1
            stage = res / (s - pole)
            stages += stage
            new_w -= stage
            if show:
                quadrant_plot(s.imag, new_w, 'After step ' + str(count))
    return new_w, stages

def nan_to_zero(a):
    r = a.copy()
    r[isnan(r)] = 0
    return r


if __name__ == '__main__':
    brune_model_stages = [(.03, .1, .2, .35, -(.2*.1)/(.2+.1)), (.01, .2, .6, .05, -(.6*.2)/(.6+.2))]
    s = 1j*linspace(0.0, 10, 501)[1:]
    # poles = [(10, 3, .3), (20, 6, .1), (30, 7, 1)]
    # z = z_init = par(*[r + zc*(s/w0 + w0/s) for zc, w0, r in poles])
    z = z_init = brune_estimate(s, brune_model_stages, 10)
    y = 1 / z
    w = imag(s)

    stages = []
    n_stages = 2
    for i in range(1, n_stages+1):
        print "\nStage", i
        print "=============="
        mask = bitwise_not(isnan(z))
        print '\nZ Preamble'
        print '---------\n'
        figure(3*i-2)
        zf, z_stages = foster_preamble(z[mask], s[mask], 4, show=True)
        z[mask] = zf
        figure(3*i-1)
        print '\nY Preamble'
        print '---------\n'
        y = 1 / z
        yf, y_stages = foster_preamble(y[mask], s[mask], 4, show=True)
        y[mask] = yf
        z = 1 / y
        figure(3*i)
        print '\nBrune Step'
        print '---------\n'
        prev_z = z
        stage, z, rf = brune_step(s, z)
        stage += z_stages, y_stages, mask
        stage_str = "r1:%s l1:%s l2:%s c2:%s l3:%s rf:%s" % (stage[:5] + (rf,))
        print stage_str
        stages.append(stage)
        quadrant_plot(w, brune_estimate(s, [stage], rf), "Z est")
        quadrant_plot(w, prev_z, "Z input")
        quadrant_plot(w, z, "Z next")
        subplot(2,2,1)
        yscale("log")
        subplot(2,2,3)
        yscale("log")
        figure(3*n_stages+1)
        quadrant_plot(w, brune_estimate(s, stages, rf), "Z stages to " + str(i))

    mask = bitwise_not(isnan(z))
    zf, _ = foster_preamble(z[mask], s[mask], 2, threshold=1)
    rf = min(zf.real)
    print 'final rf', rf

    quadrant_plot(w, brune_estimate(s, stages, rf), "Z final rf")
    quadrant_plot(w, z_init, "Z init")
    subplot(2,2,1)
    yscale("log")
    subplot(2,2,3)
    yscale("log")
    figure(3*n_stages+1)

    print [s[:5] for s in stages]
    print brune_model_stages
    show()


# def brune(y, s, order):
#     n_poles = 2*(order)
#     poles, residues, d, h = vectfit_auto(y, s, n_poles=order, n_iter=10, show=True)
#     denom_coefs = zeros(n_poles+1, dtype=complex64)
#     num_coefs = zeros(n_poles+1, dtype=complex64)
#     for n in range(n_poles+1):
#         denom_coefs[n] = sum(product(c) for c in combinations(poles, n))
#         num_coefs[n] = d * denom_coefs[n]
#         for n_prime, res in enumerate(residues):
#             poles_prime = poles[poles != poles[n_prime]]
#             num_coefs[n] += res * sum(product(c) for c in combinations(poles_prime, n))
#     denom_coefs = real(denom_coefs)
#     num_coefs = real(num_coefs)
#
#     print 'denominator', denom_coefs
#     print 'numerator', num_coefs
#
#     stages = []
#     while n_poles > 0:
#         if denom_coefs[0] == 0:
#             tmp_denom_coefs = zeros_like(denom_coefs)
#             tmp_denom_coefs[:-1] = denom_coefs[1:]
#             q, r = polydiv(num_coefs, tmp_denom_coefs)
#             assert q.shape == (1,)
#             assert r.shape == (n_poles-1)
#             stages.append("parallel", "capacitor", q[0])
#             num_coefs = r
#             denom_coefs = denom_coefs[1:]
#             n_poles -= 1
#
#         if denom_coefs[n_poles] == 0:
#             q, r = polydiv(num_coefs[::-1], denom_coefs[::-1])
#             stages.append("parallel", "inductor", q[0])
#             num_coefs = r[::-1]
#             denom_coefs = denom_coefs[1:]
#             n_poles -= 1
