"""The rows of the kernel density estimate.

The generator writes test/suite/kde_vectors.goose. The numbers come from
SciPy 1.18.1. Needs Python 3 with SciPy.
"""
from vecgen import G, f, arr

SCRIPT = "make_kde_vectors.py"


def kde_group():
    """The rows of the kernel density estimate."""
    from scipy.stats import gaussian_kde
    g = G("kde", SCRIPT)
    d1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    p1 = [0.5, 2.5, 4.3, 7.5]
    g.decl("d1", d1)
    g.decl("p1", p1)
    k = gaussian_kde(d1)
    g.raw("var k1, ck1 = statistics::kde_fit(d1[..], 1, [], statistics::KDE_SCOTT, 0.0);")
    g.raw('chk_val(bad, "kde scott factor", k1.factor, ck1, %s, 1e-12);' % f(k.factor))
    g.raw('chk_val(bad, "kde scott neff", k1.neff, ck1, %s, 1e-12);' % f(k.neff))
    g.raw("var o1: f64[>..] = [];")
    g.raw("let cp1 = statistics::kde_pdf(k1, o1, p1[..]);")
    g.raw("let w1: f64[%d] = %s;" % (len(p1), arr(k(p1))))
    g.raw('chk_list(bad, "kde pdf", o1[..], cp1, w1[..], 1e-10);')
    g.raw("var o2: f64[>..] = [];")
    g.raw("let cl1 = statistics::kde_logpdf(k1, o2, p1[..]);")
    g.raw("let w2: f64[%d] = %s;" % (len(p1), arr(k.logpdf(p1))))
    g.raw('chk_list(bad, "kde logpdf", o2[..], cl1, w2[..], 1e-10);')
    g.raw("let ib, cib = statistics::kde_integrate_box_1d(k1, 2.0, 5.0);")
    g.raw('chk_val(bad, "kde box", ib, cib, %s, 1e-10);' % f(k.integrate_box_1d(2.0, 5.0)))
    g.raw("var cv: f64[>..] = [];")
    g.raw("let ccv = statistics::kde_covariance(k1, cv);")
    g.raw("let w3: f64[1] = %s;" % arr([k.covariance[0][0]]))
    g.raw('chk_list(bad, "kde covariance", cv[..], ccv, w3[..], 1e-12);')
    # the factor of Silverman
    ks = gaussian_kde(d1, bw_method="silverman")
    g.raw("var k2, ck2 = statistics::kde_fit(d1[..], 1, [], statistics::KDE_SILVERMAN, 0.0);")
    g.raw('chk_val(bad, "kde silverman factor", k2.factor, ck2, %s, 1e-12);' % f(ks.factor))
    # the factor that the caller gives
    g.raw("var k3, ck3 = statistics::kde_fit(d1[..], 1, [], statistics::KDE_SCALAR, 0.5);")
    g.raw('chk_val(bad, "kde scalar factor", k3.factor, ck3, 0.5, 1e-12);')
    # the weighted estimate
    wts = [1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 1.0]
    g.decl("wt", wts)
    kw = gaussian_kde(d1, weights=wts)
    g.raw("var k4, ck4 = statistics::kde_fit(d1[..], 1, wt[..], statistics::KDE_SCOTT, 0.0);")
    g.raw('chk_val(bad, "kde weighted factor", k4.factor, ck4, %s, 1e-12);' % f(kw.factor))
    g.raw('chk_val(bad, "kde weighted neff", k4.neff, ck4, %s, 1e-12);' % f(kw.neff))
    g.raw("var o4: f64[>..] = [];")
    g.raw("let cp4 = statistics::kde_pdf(k4, o4, p1[..]);")
    g.raw("let w4: f64[%d] = %s;" % (len(p1), arr(kw(p1))))
    g.raw('chk_list(bad, "kde weighted pdf", o4[..], cp4, w4[..], 1e-10);')
    # a sample of two dimensions
    d2 = [0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 3.0, 1.0, 4.0, 2.0, 5.0, 0.0]
    p2 = [1.5, 1.0, 0.0, 0.0]
    g.decl("d2", d2)
    g.decl("p2", p2)
    k2d = gaussian_kde([d2[0::2], d2[1::2]])
    g.raw("var k5, ck5 = statistics::kde_fit(d2[..], 2, [], statistics::KDE_SCOTT, 0.0);")
    g.raw('chk_val(bad, "kde two factor", k5.factor, ck5, %s, 1e-12);' % f(k2d.factor))
    g.raw("var o5: f64[>..] = [];")
    g.raw("let cp5 = statistics::kde_pdf(k5, o5, p2[..]);")
    g.raw("let w5: f64[2] = %s;" % arr(k2d([[p2[0], p2[2]], [p2[1], p2[3]]])))
    g.raw('chk_list(bad, "kde two pdf", o5[..], cp5, w5[..], 1e-10);')
    g.write()


if __name__ == "__main__":
    kde_group()
