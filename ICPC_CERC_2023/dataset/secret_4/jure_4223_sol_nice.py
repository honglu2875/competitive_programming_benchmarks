import sys
from dataclasses import dataclass
from collections import Counter

@dataclass
class Entry:
    start: int
    end: int
    user: str

def parse(line):
    a, b, u = line.strip().split()
    return Entry(int(a), int(b), u)

def count(s):
    c = Counter()
    for x in s:
        c[x.user] += x.end - x.start
    return c

def diff(s1, s2):
   c = count(s2)
   c.subtract(count(s1))
   return c

s1 = []
s2 = []

s = s1
for line in sys.stdin:
    if line[0] == '-':
        s = s2
    elif line[0] == '=':
        break
    else:
        s.append(parse(line))

d = diff(s1, s2)
out = [f'{u} {d[u]:+}' for u in sorted(d) if d[u] != 0]
if out:
    print('\n'.join(out))
else:
    print("No differences found.")

