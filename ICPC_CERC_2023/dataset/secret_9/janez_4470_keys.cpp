#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cassert>
#include <random>
#include <algorithm>
#include <cstring>
#include <stack>
#include <utility>
#include <ctime>
#include <unordered_set>
#include <cmath>
using namespace std;

#define Assert assert
//#define Assert(x) 

int main()
{
    //freopen("keys97.in", "rt", stdin);
    constexpr int MaxN = 100'000, MaxM = 100'000;
    // Read the input data.
    int n, m; int ok_ = scanf("%d %d", &n, &m); Assert(ok_ == 2);
    Assert(2 <= n); Assert(n <= MaxN); Assert(2 <= m); Assert(m <= MaxM);
    typedef int RoomNo, KeyNo, CompNo;
    constexpr RoomNo bedroom = 0, outside = 1;
    struct Neigh { RoomNo u; KeyNo k; };
    vector<vector<Neigh>> neigh(n);
    for (KeyNo k = 0; k < m; ++k) {
        RoomNo a, b; ok_ = scanf("%d %d", &a, &b); Assert(ok_ == 2);
        Assert(0 <= a); Assert(a < n); Assert(0 <= b); Assert(b < n); Assert(a != b);
        neigh[a].push_back({b, k}); neigh[b].push_back({a, k}); }
    auto HasParallelEdges = [&neigh] () { 
        for (auto &N : neigh) {
            sort(N.begin(), N.end(), [] (const auto &x, const auto &y) { return x.u < y.u; });
            for (int i = 1; i < N.size(); ++i) { Assert(N[i - 1].u != N[i].u); if (N[i - 1].u == N[i].u) return true; } }
        return false; };
    Assert(! HasParallelEdges());
    // Explore the connected components if we pretend that vertex 1 doesn't exist.
    // We'll start the exploration of each component from some neighbour u of 1 and do a
    // breadth-first search; thus we get a tree of paths from u to all other vertices
    // in the component, and pred[v] is the predecessor of v on the path from u to v.
    vector<CompNo> comp(n, -1); int nComps = 0;
    vector<Neigh> pred(n, {-1, -1}); vector<RoomNo> compStarts;
    for (const Neigh &uEdge : neigh[outside]) {
        RoomNo u = uEdge.u; if (comp[u] >= 0) continue;
        comp[u] = nComps++; compStarts.emplace_back(u);
        vector<RoomNo> toDo; toDo.emplace_back(u); pred[u] = {outside, uEdge.k};
        while (! toDo.empty()) {
            RoomNo v = toDo.back(); toDo.pop_back(); Assert(comp[v] == comp[u]);
            for (const Neigh &vEdge : neigh[v]) {
                RoomNo w = vEdge.u; if (w == outside) continue;
                if (comp[w] >= 0) { Assert(comp[w] == comp[u]); continue; }
                comp[w] = comp[v]; pred[w] = {v, vEdge.k}; toDo.emplace_back(w); } } }
    // See if 1 has at least two neighbours a and b that both belong to the component which includes 0.
    for (RoomNo u = 0; u < n; ++u) if (u != outside) { Assert(comp[u] >= 0); Assert(pred[u].u >= 0); Assert(pred[u].k >= 0); }
    CompNo comp0 = comp[bedroom]; RoomNo a = compStarts[comp0], b = -1; KeyNo bKey = -1;
    for (const Neigh &e : neigh[outside]) 
        if (RoomNo u = e.u; u != a && comp[u] == comp0) { b = u; bKey = e.k; break; }
    // If yes, we have the following scenario: since 'a' is 'compStarts[comp0]',
    // the 'pred' array gives us a path from 0 to 'a' and also a path from 'b' to 'a'.
    // Let 'z' be the node closest to 'b' on this latter path that does not also lie on
    // the path from 0 to 'a'.  (Possibly 'z' may be equal to 0 or 'b', which is also fine.)
    // Then Alice can take the keys 0...z...a-1, walk 0...z, drop keys for 0...z in z,
    // then walk z...a-1; Bob can take the keys z...b-1, walk 1-b...z, pick up keys for
    // 0...z in z, then walk z...0.
    auto GetPath = [outside, &pred, &compStarts, &comp] (RoomNo u, vector<Neigh> &path) {
        path.clear(); Assert(u != outside); CompNo c = comp[u]; RoomNo u0 = u;
        while (u != outside) { path.push_back(pred[u]); u = path.back().u; } 
        Assert(path.back().u == outside);
        if (u0 == compStarts[c]) Assert(path.size() == 1);
        else Assert(path[path.size() - 2].u == compStarts[c]); };
    auto PrintKeys = [] (const vector<Neigh> & path, int iFrom = 0, int iTo = -1, bool first = true) { 
        for (int i = iFrom; i < (iTo < 0 ? path.size() : iTo); ++i) {
            if (first) first = false; else printf(" ");  printf("%d", path[i].k); } };
    auto Move = [] (RoomNo u) { printf("MOVE %d\n", u); };
    auto Done = [] () { printf("DONE\n"); };
    auto Grab = [] () { printf("GRAB\n"); };
    vector<Neigh> path0, pathB; GetPath(bedroom, path0);
    vector<KeyNo> keys;
    if (b >= 0)
    {
        GetPath(b, pathB);
        vector<int> onPath0(n, -1); for (int i = 0; i <= path0.size(); ++i) onPath0[i == 0 ? bedroom : path0[i - 1].u] = i;
        RoomNo z = -1; int zi = -1; for (int i = 0; i <= pathB.size(); ++i) if (RoomNo u = (i == 0) ? b : pathB[i - 1].u; onPath0[u] >= 0) { z = u; zi = i; break; }
        Assert(z >= 0);
        PrintKeys(path0); printf("\n"); // Alice's keys
        printf("%d", bKey); if (zi > 0) PrintKeys(pathB, 0, zi, false); printf("\n"); // Bob's keys
        for (auto &e : path0) { // Alice's moves, 0...z...a-1.
            Move(e.u); keys.emplace_back(e.k);
            if (e.u != z) continue;
            // In room 'z', Alice drops the keys for the path 0...z.
            if (! keys.empty()) { printf("DROP"); for (KeyNo k : keys) printf(" %d", k); printf("\n"); } }
        Done();
        Move(b); // Bob steps from 1 to b.
        if (b != z) for (auto &e : pathB) { // Bob's moves, b...z.
            Move(e.u); if (e.u == z) break; }
        if (z != bedroom) {
            Grab(); 
            for (int i = onPath0[z] - 1; i >= 0; --i)  // Bob's moves, z...0.
                Move(i == 0 ? bedroom : path0[i - 1].u); }
        Done();
        return 0;
    }
    // Otherwise we know that 0 can be reached from 1 through only one room, namely 'a'.
    // All other neighbours of 1 belong to other components.  Try to find two neighbours 'b' and 'c'
    // that belong to the same component.
    RoomNo c = -1; for (const Neigh &e : neigh[outside]) {
        CompNo cm = comp[e.u]; if (e.u == compStarts[cm]) continue;
        c = compStarts[cm]; b = e.u; bKey = e.k; break; }
    if (c < 0) {
        // If all neighbours of 1 belong to separate components, the problem can't be solved.
        // Alice needs the key for a-1 to get out; therefore Bob can't take it with him; 
        // therefore his path can't begin with 1-a, but with 1-b for some other neighbour b;
        // therefore he must take the key for 1-b, therefore Alice can't enter the component
        // to which b belongs and drop keys there; therefore starting his path with 1-b doesn't
        // bring Bob any closer to being able to enter room a (and from there eventually reach 0)
        // since he won't get any new keys there.
        printf("No solution\n"); return 0; }
    // Now we have the following scenario: Alice takes the keys for 0...a-1-b, walks this
    // path, drops the keys for 0...a-1 in room b, and moves back to 1.  Bob takes the keys for
    // 1-c...b, walks this path, picks up the keys for 0...a-1 in room b, then walks b...c-1-a...0.
    GetPath(b, pathB);
    printf("%d", bKey); PrintKeys(path0, 0, -1, false); printf("\n"); // Alice's keys
    PrintKeys(pathB); printf("\n"); // Bob's keys
    // Alice's moves.
    for (auto &e : path0) { Move(e.u); keys.emplace_back(e.k); } // from 0 via a to 1
    Move(b); Assert(! keys.empty()); // step into b and drop keys
    printf("DROP"); for (KeyNo k : keys) printf(" %d", k); printf("\n");
    Move(outside); Done(); // step into 1 and done
    // Bob's moves.
    for (int i = pathB.size() - 1; i >= 0; --i) Move(i == 0 ? b : pathB[i - 1].u); // from 1 via c to b
    Grab(); // pick up keys in b
    for (auto &e : pathB) Move(e.u); // from b via c to 1
    for (int i = path0.size() - 1; i >= 0; --i) Move(i == 0 ? bedroom : path0[i - 1].u); // from 1 via a to 0
    Done();
    return 0;
}