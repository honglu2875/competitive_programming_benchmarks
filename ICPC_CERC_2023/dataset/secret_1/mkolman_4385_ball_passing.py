import sys

n = int(input())
g = input()
coor = [complex(*map(int, input().split())) for _ in range(n)]

def calc(x):
    gcoor = [coor[i] for i in range(n) if g[i] == x]
    return sum(abs(a-b) for a, b in zip(gcoor, gcoor[len(gcoor)//2:]))

print(f"{calc('B') + calc('G'):.8f}")

assert 2 <= n <= 50
assert len(g) == n
assert len(list(sys.stdin)) == 0
assert all(abs(c.real) < 10000 and abs(c.imag) < 10000 for c in coor)