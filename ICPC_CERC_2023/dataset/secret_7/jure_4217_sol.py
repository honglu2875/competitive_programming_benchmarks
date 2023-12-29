import sys
 
SPEC = [2, 4, 10]
# critical lengths are 3 (2^2-1) and 18 (3 + 2^4-1)
# for 3, we can make at most (600/3) lines with 6 bits, total 1200
# for 18, we can make at most (600/18) lines with 16 bits, total 528
 
def varint_coding_factory(spec):
    size = [(1 << i)-1 for i in spec]
 
    def varint_encode(lst):
        r = []
        for n in lst:
            e = ''
            for i, s in enumerate(size):
                if n < s:
                    e += '{:0>{w}b}'.format(n, w=spec[i])
                    break
                else:
                    n -= s
                    e += '1'*spec[i]
            else:
                raise ValueError(f"{n} too large")
            r.append(e)
        return ''.join(r)
 
    def varint_decode(enc):
        r = []
        i = 0
        while i < len(enc):
            n = 0
            for j, s in enumerate(spec):
                p = int(enc[i:i+s], 2)
                i += s
                if p < size[j]:
                    r.append(n + p)
                    break
                else:
                    n += p
            else:
                raise ValueError(f"Cannot parse {enc}")
        return r
 
    return varint_encode, varint_decode
 
encode, decode = varint_coding_factory(SPEC)
 
cmd = sys.stdin.readline().strip()
if cmd == 'ENCODE':
    e = {}
    for line in sys.stdin:
        r, c = line.strip().split(": ")
        cs = c.split()
        e[r] = cs
 
 
    lvl = []
    pd = -1
    def pre(r, d=0):
        global pd
        print(r)
        lvl.append(pd-d+1)
        pd = d
        for c in e.get(r, []):
            pre(c, d+1)
 
    for root in e:
        pre(root)
        break
 
    enc = encode(lvl)
    print(enc)
 
elif cmd == 'DECODE':
    employees = [line.strip() for line in sys.stdin]
    ds = decode(employees.pop())
    root = employees[0]
    s = {root: []}
    st = [root]
    for e, d in zip(employees[1:], ds[1:]):
        for i in range(d): st.pop()
        s[e] = []
        s[st[-1]].append(e)
        st.append(e)
 
    for e, l in list(s.items()):
        if len(l) == 0:
            del s[e]
 
    for r, cs in s.items():
        print(f'{r}: {" ".join(cs)}')
 
else:

    assert False, cmd
