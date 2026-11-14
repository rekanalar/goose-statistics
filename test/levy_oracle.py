"""The oracle of the stable law, for the vector file of the law.

The oracle gives the density and the cumulative function of the law by
the Fourier representation of the law. The representation is a different
formula from the integral of Nolan that the package uses. The two agree
to about 1e-15 at every point of the grid.

A shape of one has its own form of the representation. The tangent of
the general form goes to infinity at a shape of one, so the general form
does not reach it. The limit of the characteristic function gives the
form below.

Needs: Python 3 and mpmath. Run it from the package directory with the
system python3.
"""
import mpmath as mp

mp.mp.dps = 30

# The panels of the integral of the Fourier form. The integrand falls
# off like the exponential of the power of the angle, so a geometric set
# of panels covers the range.
CUTS = [mp.mpf(2) ** k for k in range(-4, 12)]

# The panels of the integral for a shape of one. The integrand falls off
# like the exponential of the angle.
CUTS_ONE = [mp.mpf(2) ** k for k in range(0, 9)]


def characteristic(t, alpha, beta):
    """The characteristic function of the law at a positive angle.

    The function of a shape of one is the limit of the general form.
    """
    if alpha == 1:
        return mp.e ** (-t) * mp.e ** (-1j * beta * 2 / mp.pi * t * mp.log(t))
    tb = mp.tan(mp.pi * alpha / 2)
    return mp.e ** (-t ** alpha * (1 - 1j * beta * tb))


def pdf(x, alpha, beta):
    """The density of the law, by the Fourier representation."""
    x, alpha, beta = mp.mpf(x), mp.mpf(alpha), mp.mpf(beta)
    cuts = CUTS_ONE if alpha == 1 else CUTS

    def integrand(t):
        return mp.re(characteristic(t, alpha, beta) * mp.e ** (-1j * x * t))

    return mp.quad(integrand, [0] + cuts) / mp.pi


def cdf(x, alpha, beta):
    """The cumulative function of the law, by the Fourier representation.

    The integrand carries the factor 1 over the imaginary unit t, so
    the integral stays at the lower end of the range.
    """
    x, alpha, beta = mp.mpf(x), mp.mpf(alpha), mp.mpf(beta)
    cuts = CUTS_ONE if alpha == 1 else CUTS

    def integrand(t):
        return mp.re(characteristic(t, alpha, beta) * mp.e ** (-1j * x * t)
                     / (-1j * t))

    return mp.mpf("0.5") + mp.quad(integrand, [0] + cuts) / mp.pi


def pair(job):
    """The density and the cumulative function of one point."""
    a, b, x = job
    return (float(pdf(x, a, b)), float(cdf(x, a, b)))
