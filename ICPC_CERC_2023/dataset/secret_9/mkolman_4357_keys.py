import sys
from collections import defaultdict as dd

sys.setrecursionlimit(10**6)

def find_path(graph, visited=None, current=0, end=1, parent=None):
    if visited is None:
        visited = set()
    visited.add(current)
    for n in graph[current]:
        if n == parent:
            continue
        if n == end:
            return [end, current]
        if n not in visited:
            r = find_path(graph, visited, n, end, current)
            if r is not None:
                r.append(current)
                return r
    return None

def move(path, extra_at, extra):
    for room in path[1:]:
        print(f"MOVE {room}")
        if extra_at == room:
            print(extra)
    print("DONE")

def cycle_overlap_at(path, cycle):
    return next(filter(set(path[:-1]).__contains__, cycle), None)

def keys(graph, path):
    return ' '.join(set(graph[i][j] for i, j in zip(path, path[1:])))

def main():
    n, m = map(int, input().split())
    assert 2 <= n <= 10**5 
    assert 2 <= m <= 10**5 
    graph = dd(dict)
    for i in range(m):
        a, b = map(int, input().split())
        assert 0 <= a < n
        assert 0 <= b < n
        assert a != b
        assert b not in graph[a]
        assert a not in graph[b]
        graph[a][b] = str(i)
        graph[b][a] = str(i)
    path = find_path(graph)
    if path is None:
        print("No solution")
        return
    cycle = find_path(graph, current=1, parent=path[1])
    if cycle is None:
        print("No solution")
        return
    path = path[::-1]
    cycle = cycle[::-1]
    
    overlap_at = cycle_overlap_at(path, cycle)

    if overlap_at is None:
        print(keys(graph, path+[cycle[-2]]))
        print(keys(graph, cycle[:-1]))
        move(path+cycle[-2:], cycle[-2], f"DROP {keys(graph, path)}")
        move(cycle[:-1] + cycle[-3::-1] + path[-2::-1], cycle[-2], "GRAB")
    else:
        path_idx = path.index(overlap_at)
        cycle_idx = cycle.index(overlap_at)
        print(keys(graph, path))
        print(keys(graph, cycle[:cycle_idx+1]))
        move(path, overlap_at, f"DROP {' '.join(graph[i][j] for i, j in zip(path[:path_idx], path[1:]))}")
        move(cycle[:cycle_idx] + path[path_idx::-1], overlap_at, "GRAB")
    
if __name__ == '__main__':
    main()
