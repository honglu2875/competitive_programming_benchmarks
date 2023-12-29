import sys
from collections import defaultdict as dd
from string import ascii_lowercase

oncall = [dd(int)]
prev = '0'
while True:
    line = input()
    if line == "------":
        oncall.append(dd(int))
        prev = '0'
        continue
    if line == "======":
        break
    s, e, n = line.split(" ")
    assert s.strip() == s and e.strip() == e and n.strip() == n
    assert int(s) < int(e)
    assert int(e) <= 1000
    assert 3 <= len(n) <= 20
    assert all(c in ascii_lowercase for c in n)
    assert s == prev
    oncall[-1][n] += int(e) - int(s)
    prev = e

assert len(oncall) == 2
assert len(list(sys.stdin)) == 0
change = False
for name in sorted(set(list(oncall[0].keys())+list(oncall[1].keys()))):
    diff = oncall[1][name] - oncall[0][name]
    if diff:
        change = True
        print(f"{name} {diff:+d}")
if not change:
    print("No differences found.")