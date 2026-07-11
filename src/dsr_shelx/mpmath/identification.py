"""
Implements the PSLQ algorithm for integer relation detection,
and derivative algorithms for constant recognition.
"""

from .libmp.backend import xrange
from .libmp import int_types, sqrt_fixed

# round to nearest integer (can be done more elegantly...)
def round_fixed(x, prec):
    return ((x + (1<<(prec-1))) >> prec) << prec

class IdentificationMethods(object):
    pass


def pslq(ctx, x, tol=None, maxcoeff=1000, maxsteps=100, verbose=False):
    r"""
    Given a vector of real numbers `x = [x_0, x_1, ..., x_n]`, ``pslq(x)``
    uses the PSLQ algorithm to find a list of integers
    `[c_0, c_1, ..., c_n]` such that

    .. math ::

        |c_1 x_1 + c_2 x_2 + ... + c_n x_n| < \mathrm{tol}

    and such that `\max |c_k| < \mathrm{maxcoeff}`. If no such vector
    exists, :func:`~mpmath.pslq` returns ``None``. The tolerance defaults to
    3/4 of the working precision.

    **Examples**

    Find rational approximations for `\pi`::

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> pslq([-1, pi], tol=0.01)
        [22, 7]
        >>> pslq([-1, pi], tol=0.001)
        [355, 113]
        >>> mpf(22)/7; mpf(355)/113; +pi
        3.14285714285714
        3.14159292035398
        3.14159265358979

    Pi is not a rational number with denominator less than 1000::

        >>> pslq([-1, pi])
        >>>

    To within the standard precision, it can however be approximated
    by at least one rational number with denominator less than `10^{12}`::

        >>> p, q = pslq([-1, pi], maxcoeff=10**12)
        >>> print(p); print(q)
        238410049439
        75888275702
        >>> mpf(p)/q
        3.14159265358979

    The PSLQ algorithm can be applied to long vectors. For example,
    we can investigate the rational (in)dependence of integer square
    roots::

        >>> mp.dps = 30
        >>> pslq([sqrt(n) for n in range(2, 5+1)])
        >>>
        >>> pslq([sqrt(n) for n in range(2, 6+1)])
        >>>
        >>> pslq([sqrt(n) for n in range(2, 8+1)])
        [2, 0, 0, 0, 0, 0, -1]

    **Machin formulas**

    A famous formula for `\pi` is Machin's,

    .. math ::

        \frac{\pi}{4} = 4 \operatorname{acot} 5 - \operatorname{acot} 239

    There are actually infinitely many formulas of this type. Two
    others are

    .. math ::

        \frac{\pi}{4} = \operatorname{acot} 1

        \frac{\pi}{4} = 12 \operatorname{acot} 49 + 32 \operatorname{acot} 57
            + 5 \operatorname{acot} 239 + 12 \operatorname{acot} 110443

    We can easily verify the formulas using the PSLQ algorithm::

        >>> mp.dps = 30
        >>> pslq([pi/4, acot(1)])
        [1, -1]
        >>> pslq([pi/4, acot(5), acot(239)])
        [1, -4, 1]
        >>> pslq([pi/4, acot(49), acot(57), acot(239), acot(110443)])
        [1, -12, -32, 5, -12]

    We could try to generate a custom Machin-like formula by running
    the PSLQ algorithm with a few inverse cotangent values, for example
    acot(2), acot(3) ... acot(10). Unfortunately, there is a linear
    dependence among these values, resulting in only that dependence
    being detected, with a zero coefficient for `\pi`::

        >>> pslq([pi] + [acot(n) for n in range(2,11)])
        [0, 1, -1, 0, 0, 0, -1, 0, 0, 0]

    We get better luck by removing linearly dependent terms::

        >>> pslq([pi] + [acot(n) for n in range(2,11) if n not in (3, 5)])
        [1, -8, 0, 0, 4, 0, 0, 0]

    In other words, we found the following formula::

        >>> 8*acot(2) - 4*acot(7)
        3.14159265358979323846264338328
        >>> +pi
        3.14159265358979323846264338328

    **Algorithm**

    This is a fairly direct translation to Python of the pseudocode given by
    David Bailey, "The PSLQ Integer Relation Algorithm":
    http://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html

    The present implementation uses fixed-point instead of floating-point
    arithmetic, since this is significantly (about 7x) faster.
    """

    n = len(x)
    if n < 2:
        raise ValueError("n cannot be less than 2")

    # At too low precision, the algorithm becomes meaningless
    prec = ctx.prec
    if prec < 53:
        raise ValueError("prec cannot be less than 53")

    if verbose and prec // max(2,n) < 5:
        print("Warning: precision for PSLQ may be too low")

    target = int(prec * 0.75)

    if tol is None:
        tol = ctx.mpf(2)**(-target)
    else:
        tol = ctx.convert(tol)

    extra = 60
    prec += extra

    if verbose:
        print("PSLQ using prec %i and tol %s" % (prec, ctx.nstr(tol)))

    tol = ctx.to_fixed(tol, prec)
    assert tol

    # Convert to fixed-point numbers. The dummy None is added so we can
    # use 1-based indexing. (This just allows us to be consistent with
    # Bailey's indexing. The algorithm is 100 lines long, so debugging
    # a single wrong index can be painful.)
    x = [None] + [ctx.to_fixed(ctx.mpf(xk), prec) for xk in x]

    # Sanity check on magnitudes
    minx = min(abs(xx) for xx in x[1:])
    if not minx:
        raise ValueError("PSLQ requires a vector of nonzero numbers")
    if minx < tol//100:
        if verbose:
            print("STOPPING: (one number is too small)")
        return None

    g = sqrt_fixed((4<<prec)//3, prec)
    A = {}
    B = {}
    H = {}
    # Initialization
    # step 1
    for i in xrange(1, n+1):
        for j in xrange(1, n+1):
            A[i,j] = B[i,j] = (i==j) << prec
            H[i,j] = 0
    # step 2
    s = [None] + [0] * n
    for k in xrange(1, n+1):
        t = 0
        for j in xrange(k, n+1):
            t += (x[j]**2 >> prec)
        s[k] = sqrt_fixed(t, prec)
    t = s[1]
    y = x[:]
    for k in xrange(1, n+1):
        y[k] = (x[k] << prec) // t
        s[k] = (s[k] << prec) // t
    # step 3
    for i in xrange(1, n+1):
        for j in xrange(i+1, n):
            H[i,j] = 0
        if i <= n-1:
            if s[i]:
                H[i,i] = (s[i+1] << prec) // s[i]
            else:
                H[i,i] = 0
        for j in range(1, i):
            sjj1 = s[j]*s[j+1]
            if sjj1:
                H[i,j] = ((-y[i]*y[j])<<prec)//sjj1
            else:
                H[i,j] = 0
    # step 4
    for i in xrange(2, n+1):
        for j in xrange(i-1, 0, -1):
            #t = floor(H[i,j]/H[j,j] + 0.5)
            if H[j,j]:
                t = round_fixed((H[i,j] << prec)//H[j,j], prec)
            else:
                #t = 0
                continue
            y[j] = y[j] + (t*y[i] >> prec)
            for k in xrange(1, j+1):
                H[i,k] = H[i,k] - (t*H[j,k] >> prec)
            for k in xrange(1, n+1):
                A[i,k] = A[i,k] - (t*A[j,k] >> prec)
                B[k,j] = B[k,j] + (t*B[k,i] >> prec)
    # Main algorithm
    for REP in range(maxsteps):
        # Step 1
        m = -1
        szmax = -1
        for i in range(1, n):
            h = H[i,i]
            sz = (g**i * abs(h)) >> (prec*(i-1))
            if sz > szmax:
                m = i
                szmax = sz
        # Step 2
        y[m], y[m+1] = y[m+1], y[m]
        for i in xrange(1,n+1): H[m,i], H[m+1,i] = H[m+1,i], H[m,i]
        for i in xrange(1,n+1): A[m,i], A[m+1,i] = A[m+1,i], A[m,i]
        for i in xrange(1,n+1): B[i,m], B[i,m+1] = B[i,m+1], B[i,m]
        # Step 3
        if m <= n - 2:
            t0 = sqrt_fixed((H[m,m]**2 + H[m,m+1]**2)>>prec, prec)
            # A zero element probably indicates that the precision has
            # been exhausted. XXX: this could be spurious, due to
            # using fixed-point arithmetic
            if not t0:
                break
            t1 = (H[m,m] << prec) // t0
            t2 = (H[m,m+1] << prec) // t0
            for i in xrange(m, n+1):
                t3 = H[i,m]
                t4 = H[i,m+1]
                H[i,m] = (t1*t3+t2*t4) >> prec
                H[i,m+1] = (-t2*t3+t1*t4) >> prec
        # Step 4
        for i in xrange(m+1, n+1):
            for j in xrange(min(i-1, m+1), 0, -1):
                try:
                    t = round_fixed((H[i,j] << prec)//H[j,j], prec)
                # Precision probably exhausted
                except ZeroDivisionError:
                    break
                y[j] = y[j] + ((t*y[i]) >> prec)
                for k in xrange(1, j+1):
                    H[i,k] = H[i,k] - (t*H[j,k] >> prec)
                for k in xrange(1, n+1):
                    A[i,k] = A[i,k] - (t*A[j,k] >> prec)
                    B[k,j] = B[k,j] + (t*B[k,i] >> prec)
        # Until a relation is found, the error typically decreases
        # slowly (e.g. a factor 1-10) with each step TODO: we could
        # compare err from two successive iterations. If there is a
        # large drop (several orders of magnitude), that indicates a
        # "high quality" relation was detected. Reporting this to
        # the user somehow might be useful.
        best_err = maxcoeff<<prec
        for i in xrange(1, n+1):
            err = abs(y[i])
            # Maybe we are done?
            if err < tol:
                # We are done if the coefficients are acceptable
                vec = [int(round_fixed(B[j,i], prec) >> prec) for j in \
                range(1,n+1)]
                if max(abs(v) for v in vec) < maxcoeff:
                    if verbose:
                        print("FOUND relation at iter %i/%i, error: %s" % \
                            (REP, maxsteps, ctx.nstr(err / ctx.mpf(2)**prec, 1)))
                    return vec
            best_err = min(err, best_err)
        # Calculate a lower bound for the norm. We could do this
        # more exactly (using the Euclidean norm) but there is probably
        # no practical benefit.
        recnorm = max(abs(h) for h in H.values())
        if recnorm:
            norm = ((1 << (2*prec)) // recnorm) >> prec
            norm //= 100
        else:
            norm = ctx.inf
        if verbose:
            print("%i/%i:  Error: %8s   Norm: %s" % \
                (REP, maxsteps, ctx.nstr(best_err / ctx.mpf(2)**prec, 1), norm))
        if norm >= maxcoeff:
            break
    if verbose:
        print("CANCELLING after step %i/%i." % (REP, maxsteps))
        print("Could not find an integer relation. Norm bound: %s" % norm)
    return None

def findpoly(ctx, x, n=1, **kwargs):
    r"""
    ``findpoly(x, n)`` returns the coefficients of an integer
    polynomial `P` of degree at most `n` such that `P(x) \approx 0`.
    If no polynomial having `x` as a root can be found,
    :func:`~mpmath.findpoly` returns ``None``.

    :func:`~mpmath.findpoly` works by successively calling :func:`~mpmath.pslq` with
    the vectors `[1, x]`, `[1, x, x^2]`, `[1, x, x^2, x^3]`, ...,
    `[1, x, x^2, .., x^n]` as input. Keyword arguments given to
    :func:`~mpmath.findpoly` are forwarded verbatim to :func:`~mpmath.pslq`. In
    particular, you can specify a tolerance for `P(x)` with ``tol``
    and a maximum permitted coefficient size with ``maxcoeff``.

    For large values of `n`, it is recommended to run :func:`~mpmath.findpoly`
    at high precision; preferably 50 digits or more.

    **Examples**

    By default (degree `n = 1`), :func:`~mpmath.findpoly` simply finds a linear
    polynomial with a rational root::

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> findpoly(0.7)
        [-10, 7]

    The generated coefficient list is valid input to ``polyval`` and
    ``polyroots``::

        >>> nprint(polyval(findpoly(phi, 2), phi), 1)
        -2.0e-16
        >>> for r in polyroots(findpoly(phi, 2)):
        ...     print(r)
        ...
        -0.618033988749895
        1.61803398874989

    Numbers of the form `m + n \sqrt p` for integers `(m, n, p)` are
    solutions to quadratic equations. As we find here, `1+\sqrt 2`
    is a root of the polynomial `x^2 - 2x - 1`::

        >>> findpoly(1+sqrt(2), 2)
        [1, -2, -1]
        >>> findroot(lambda x: x**2 - 2*x - 1, 1)
        2.4142135623731

    Despite only containing square roots, the following number results
    in a polynomial of degree 4::

        >>> findpoly(sqrt(2)+sqrt(3), 4)
        [1, 0, -10, 0, 1]

    In fact, `x^4 - 10x^2 + 1` is the *minimal polynomial* of
    `r = \sqrt 2 + \sqrt 3`, meaning that a rational polynomial of
    lower degree having `r` as a root does not exist. Given sufficient
    precision, :func:`~mpmath.findpoly` will usually find the correct
    minimal polynomial of a given algebraic number.

    **Non-algebraic numbers**

    If :func:`~mpmath.findpoly` fails to find a polynomial with given
    coefficient size and tolerance constraints, that means no such
    polynomial exists.

    We can verify that `\pi` is not an algebraic number of degree 3 with
    coefficients less than 1000::

        >>> mp.dps = 15
        >>> findpoly(pi, 3)
        >>>

    It is always possible to find an algebraic approximation of a number
    using one (or several) of the following methods:

        1. Increasing the permitted degree
        2. Allowing larger coefficients
        3. Reducing the tolerance

    One example of each method is shown below::

        >>> mp.dps = 15
        >>> findpoly(pi, 4)
        [95, -545, 863, -183, -298]
        >>> findpoly(pi, 3, maxcoeff=10000)
        [836, -1734, -2658, -457]
        >>> findpoly(pi, 3, tol=1e-7)
        [-4, 22, -29, -2]

    It is unknown whether Euler's constant is transcendental (or even
    irrational). We can use :func:`~mpmath.findpoly` to check that if is
    an algebraic number, its minimal polynomial must have degree
    at least 7 and a coefficient of magnitude at least 1000000::

        >>> mp.dps = 200
        >>> findpoly(euler, 6, maxcoeff=10**6, tol=1e-100, maxsteps=1000)
        >>>

    Note that the high precision and strict tolerance is necessary
    for such high-degree runs, since otherwise unwanted low-accuracy
    approximations will be detected. It may also be necessary to set
    maxsteps high to prevent a premature exit (before the coefficient
    bound has been reached). Running with ``verbose=True`` to get an
    idea what is happening can be useful.
    """
    x = ctx.mpf(x)
    if n < 1:
        raise ValueError("n cannot be less than 1")
    if x == 0:
        return [1, 0]
    xs = [ctx.mpf(1)]
    for i in range(1,n+1):
        xs.append(x**i)
        a = ctx.pslq(xs, **kwargs)
        if a is not None:
            return a[::-1]

def fracgcd(p, q):
    x, y = p, q
    while y:
        x, y = y, x % y
    if x != 1:
        p //= x
        q //= x
    if q == 1:
        return p
    return p, q

def pslqstring(r, constants):
    q = r[0]
    r = r[1:]
    s = []
    for i in range(len(r)):
        p = r[i]
        if p:
            z = fracgcd(-p,q)
            cs = constants[i][1]
            if cs == '1':
                cs = ''
            else:
                cs = '*' + cs
            if isinstance(z, int_types):
                if z > 0: term = str(z) + cs
                else:     term = ("(%s)" % z) + cs
            else:
                term = ("(%s/%s)" % z) + cs
            s.append(term)
    s = ' + '.join(s)
    if '+' in s or '*' in s:
        s = '(' + s + ')'
    return s or '0'

def prodstring(r, constants):
    q = r[0]
    r = r[1:]
    num = []
    den = []
    for i in range(len(r)):
        p = r[i]
        if p:
            z = fracgcd(-p,q)
            cs = constants[i][1]
            if isinstance(z, int_types):
                if abs(z) == 1: t = cs
                else:           t = '%s**%s' % (cs, abs(z))
                ([num,den][z<0]).append(t)
            else:
                t = '%s**(%s/%s)' % (cs, abs(z[0]), z[1])
                ([num,den][z[0]<0]).append(t)
    num = '*'.join(num)
    den = '*'.join(den)
    if num and den: return "(%s)/(%s)" % (num, den)
    if num: return num
    if den: return "1/(%s)" % den

def quadraticstring(ctx,t,a,b,c):
    if c < 0:
        a,b,c = -a,-b,-c
    u1 = (-b+ctx.sqrt(b**2-4*a*c))/(2*c)
    u2 = (-b-ctx.sqrt(b**2-4*a*c))/(2*c)
    if abs(u1-t) < abs(u2-t):
        if b:  s = '((%s+sqrt(%s))/%s)' % (-b,b**2-4*a*c,2*c)
        else:  s = '(sqrt(%s)/%s)' % (-4*a*c,2*c)
    else:
        if b:  s = '((%s-sqrt(%s))/%s)' % (-b,b**2-4*a*c,2*c)
        else:  s = '(-sqrt(%s)/%s)' % (-4*a*c,2*c)
    return s

# Transformation y = f(x,c), with inverse function x = f(y,c)
# The third entry indicates whether the transformation is
# redundant when c = 1
transforms = [
  (lambda ctx,x,c: x*c, '$y/$c', 0),
  (lambda ctx,x,c: x/c, '$c*$y', 1),
  (lambda ctx,x,c: c/x, '$c/$y', 0),
  (lambda ctx,x,c: (x*c)**2, 'sqrt($y)/$c', 0),
  (lambda ctx,x,c: (x/c)**2, '$c*sqrt($y)', 1),
  (lambda ctx,x,c: (c/x)**2, '$c/sqrt($y)', 0),
  (lambda ctx,x,c: c*x**2, 'sqrt($y)/sqrt($c)', 1),
  (lambda ctx,x,c: x**2/c, 'sqrt($c)*sqrt($y)', 1),
  (lambda ctx,x,c: c/x**2, 'sqrt($c)/sqrt($y)', 1),
  (lambda ctx,x,c: ctx.sqrt(x*c), '$y**2/$c', 0),
  (lambda ctx,x,c: ctx.sqrt(x/c), '$c*$y**2', 1),
  (lambda ctx,x,c: ctx.sqrt(c/x), '$c/$y**2', 0),
  (lambda ctx,x,c: c*ctx.sqrt(x), '$y**2/$c**2', 1),
  (lambda ctx,x,c: ctx.sqrt(x)/c, '$c**2*$y**2', 1),
  (lambda ctx,x,c: c/ctx.sqrt(x), '$c**2/$y**2', 1),
  (lambda ctx,x,c: ctx.exp(x*c), 'log($y)/$c', 0),
  (lambda ctx,x,c: ctx.exp(x/c), '$c*log($y)', 1),
  (lambda ctx,x,c: ctx.exp(c/x), '$c/log($y)', 0),
  (lambda ctx,x,c: c*ctx.exp(x), 'log($y/$c)', 1),
  (lambda ctx,x,c: ctx.exp(x)/c, 'log($c*$y)', 1),
  (lambda ctx,x,c: c/ctx.exp(x), 'log($c/$y)', 0),
  (lambda ctx,x,c: ctx.ln(x*c), 'exp($y)/$c', 0),
  (lambda ctx,x,c: ctx.ln(x/c), '$c*exp($y)', 1),
  (lambda ctx,x,c: ctx.ln(c/x), '$c/exp($y)', 0),
  (lambda ctx,x,c: c*ctx.ln(x), 'exp($y/$c)', 1),
  (lambda ctx,x,c: ctx.ln(x)/c, 'exp($c*$y)', 1),
  (lambda ctx,x,c: c/ctx.ln(x), 'exp($c/$y)', 0),
]

def identify(ctx, x, constants=[], tol=None, maxcoeff=1000, full=False,
    verbose=False):
    return None

IdentificationMethods.pslq = pslq
IdentificationMethods.findpoly = findpoly
IdentificationMethods.identify = identify


if __name__ == '__main__':
    import doctest
    doctest.testmod()
