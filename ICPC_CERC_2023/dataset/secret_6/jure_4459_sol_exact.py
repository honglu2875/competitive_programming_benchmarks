import math

def solve(za, zb, zc, r):
    if r == 0:
        return abs(za-zc) + abs(zb-zc)

    # one point inside, one outside
    a_in = abs(zc-za) <= r
    b_in = abs(zc-zb) <= r
    if a_in or b_in:
        return abs(za-zb)

    if za == zb:
        return abs(za-zc) + abs(zb-zc) - 2*r

    # line crosses the circle
    b = zb-za
    pc = zc-za
    proj = (pc.real*b.real+pc.imag*b.imag)/(b.real**2+b.imag**2)
    dist = abs(pc.real*b.imag - pc.imag*b.real)/abs(b)
    if 0 <= proj <= 1 and dist < r:
        return abs(za-zb)

    # General problem: reduce to unit circle.
    za -= zc; zb -= zc
    za /= r; zb /= r;

    p = za*zb
    s = za+zb

    # https://de.wikipedia.org/wiki/Alhazensches_Problem
    poly = [p.conjugate(), -s.conjugate(), 0, s, -p]
    # import numpy as np
    # opts = np.roots([p.conjugate(), -s.conjugate(), 0, s, -p])
    opts = solve_poly(poly[::-1])
    return r*min(abs(za-u) + abs(zb-u) for u in opts if abs(abs(u)-1) < 1e-6)


################ substitue for numpy.roots ##############
# adapted (corrected) from: https://github.com/yairchu/quartic

import math

def solve_poly(coefs):
    if not coefs:
        return []
    if abs(coefs[0]) == 0:
        return [0] + solve_poly(coefs[1:])
    if abs(coefs[-1]) == 0:
        return solve_poly(coefs[:-1])
    return solve_normalized_poly([x / coefs[-1] for x in coefs[:-1]])


def solve_normalized_poly(coefs):
    degree = len(coefs)
    shift = -coefs[-1] / degree
    return [
        x + shift
        for x in solve_depressed_poly(shifted_coefs(shift, coefs + [1])[: degree - 1])
    ]


def shifted_coefs(shift, coefs):
    result = [0]*len(coefs)
    for coef, bin in zip(coefs, binomials()):
        x = 1
        for i, b in enumerate(bin):
            result[len(bin) - 1 - i] += coef*b*x
            x *= shift
    return result


def binomials():
    cur = [1]
    while True:
        yield cur
        cur = [x + y for x, y in zip([0] + cur, cur + [0])]


def solve_depressed_poly(coefs):
    if not coefs:
        # Poly is: x + 0 = 0
        return [0]
    if abs(coefs[0]) == 0:
        return solve_depressed_poly(coefs[1:])
    if len(coefs) == 1:
        # Quadratic
        return sqrts(-coefs[0])
    if len(coefs) == 2:
        return solve_depressed_cubic(coefs[0], coefs[1])
    if len(coefs) == 3:
        return solve_depressed_quartic(coefs[0], coefs[1], coefs[2])
    raise ValueError("unsupported polynomial degree")


# Based on http://en.wikipedia.org/wiki/Quartic_function#Quick_and_memorable_solution_from_first_principles
def solve_depressed_quartic(e, d, c):
    if abs(d) == 0:  # biquadratic
        return [s for x in solve_poly([e, c, 1]) for s in sqrts(x)]
    sols = solve_poly([-d*d, c*c - 4*e, 2*c, 1])
    p = max(sols, key=abs) ** 0.5
    return solve_poly([c + p*p - d/p, 2*p, 2]) + solve_poly(
        [c + p*p + d/p, -2*p, 2]
    )

def sqrts(x):
    s = x**0.5
    return [-s, s]

third = 1/3.0
def cubic_root(x):
    if x.real >= 0:
        return x**third
    else:
        return -(-x)**third

# Based on http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method
third_root_unity = math.e ** (math.pi*2j / 3)
def solve_depressed_cubic(q, p):
    if abs(p) == 0:
        r = -cubic_root(q)
        return [r, r*third_root_unity, r*third_root_unity**2]
    else:
        u = cubic_root(max(solve_poly([-p*p*p / 27, q, 1]), key=abs))
        return [v - p/3/v for v in [u, u*third_root_unity, u*third_root_unity**2]]

t = int(input())
for c in range(t):
    xa, ya, xb, yb, xc, yc, r = map(int, input().split())
    print(solve(complex(xa, ya), complex(xb, yb), complex(xc, yc), r))