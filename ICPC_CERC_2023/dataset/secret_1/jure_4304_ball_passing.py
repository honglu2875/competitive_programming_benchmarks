n = int(input())
s = input().strip()
b, g = [], []
for i in range(n):
    x, y = map(float, input().split())
    [b, g][s[i] == 'G'].append(complex(x, y))

def solve(s):
    # For every convex quadrilateral it's better to take diagonals than a pair of opposing sides.
    # This can be proven by applying triangle inequality to 4 triangles created by intersecting
    # diagonals (intersection exists due to convexity).
    # For 2n-gon, any non-intersecting pair of diagonals can be improved by making it intersect (can
    # always be achieved because of convexity), so that points on opposite sides end up paired
    # together.
    k = len(s) // 2
    return sum(abs(s[i] - s[i-k]) for i in range(k))

print(solve(b) + solve(g))
