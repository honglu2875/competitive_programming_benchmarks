import sys
from collections import defaultdict as dd

sys.setrecursionlimit(10000000)

def main():
    cmd = input()
    if cmd == "ENCODE":
        encode()
    elif cmd == "DECODE":
        decode()
 
def encode():
    graph = dd(list)
    ceo = None
    for line in sys.stdin:
        boss, reports = line.split(': ')
        assert boss == boss.strip()
        graph[boss] = reports.strip().split(' ')
        if ceo is None:
            ceo = boss
    print(cypher(graph, ceo).rstrip('1'))
 
def cypher(graph, person):
     print(person)
     return ''.join(f'0{cypher(graph, employee)}1' for employee in graph[person])
 
def decode():
    graph = dd(list)
    managers = {}
    ceo, *names, code = map(str.strip, sys.stdin)
    current = ceo
    code = iter(code)
    for name in names:
        while next(code) == '1':
            current = managers[current]
        graph[current].append(name)
        managers[name] = current
        current = name
    q = [ceo]
    prevq = []
    while q:
        for employee in q:
            if reports := graph.get(employee):
                prevq += graph[employee]
                print(f'{employee}: {" ".join(reports)}')
        q, prevq = prevq, []
 
if __name__ == '__main__':
    main()
