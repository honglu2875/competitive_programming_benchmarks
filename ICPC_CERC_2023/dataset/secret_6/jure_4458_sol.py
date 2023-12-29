import math

fi = (math.sqrt(5) + 1) / 2

def gss(f, a, b, tol=1e-5):
    c = b - (b-a) / fi
    d = a + (b-a) / fi
    while abs(b-a) > tol:
        if f(c) < f(d):  # f(c) > f(d) to find the maximum
            b = d
        else:
            a = c

        # We recompute both c and d here to avoid loss of precision.
        c = b - (b-a) / fi
        d = a + (b-a) / fi
    return (a+b) / 2

def solve(a: complex, b: complex, c: complex, r: float) -> float:
    if r == 0:
        return abs(c-a) + abs(c-b)

    a_in = abs(c-a) <= r
    b_in = abs(c-b) <= r
    if a_in or b_in:
        return abs(a-b)

    if a == b:
        return abs(a-c) + abs(b-c) - 2*r

    # line crosses the circle
    o = b-a
    pc = c-a
    proj = (pc.real*o.real+pc.imag*o.imag)/(o.real**2+o.imag**2)
    dist = abs(pc.real*o.imag - pc.imag*o.real)/abs(o)
    if 0 <= proj <= 1 and dist < r:
        return abs(a-b)

    # general case
    d = b-a
    def dist(lam: float):
        t = a + lam*d
        er = t-c
        er /= abs(er)
        clam = c + r*er
        return abs(a-clam) + abs(b-clam)

    lam = gss(dist, 0, 1, 1e-8)
    return dist(lam)

t = int(input())
assert 1 <= t <= 1000
for _ in range(t):
    xa, ya, xb, yb, xc, yc, r = map(int, input().split())
    assert -1000 <= xa <= 1000
    assert -1000 <= ya <= 1000
    assert -1000 <= xb <= 1000
    assert -1000 <= yb <= 1000
    assert -1000 <= xc <= 1000
    assert -1000 <= yc <= 1000
    assert 0 <= r <= 1000
    print(solve(complex(xa, ya), complex(xb, yb), complex(xc, yc), r))

try:
    input()
    raise ValueError("More??")
except EOFError:
    pass