#include <cstdint>
#include <array>
#include <cstdio>
#include <cassert>
using namespace std;
#define Assert assert

enum { MaxI = 200, MaxC = 200, MaxT = 200, MaxPrice = 1'000'000'000, MaxAmount = 1'000'000'000 };
typedef int_fast64_t myint;

int main()
{
    //freopen("cakes99.in", "rt", stdin);
    //freopen("cakes\\test0.in", "rt", stdin);
    int I, C, T;
    int cakePrice[MaxC], ingPrice[MaxI], toolPrice[MaxT], nToolsNeeded[MaxC], nNeededBy[MaxT];
    int needLists[MaxC][MaxT], neededByLists[MaxT][MaxC];
    int ok_; ok_ = scanf("%d %d %d\n", &I, &C, &T); Assert(ok_ == 3);
    Assert(1 <= I); Assert(I <= MaxI);
    Assert(1 <= C); Assert(C <= MaxC);
    Assert(1 <= T); Assert(T <= MaxT);
    for (int ci = 0; ci < C; ++ci) { ok_ = scanf("%d", &cakePrice[ci]); Assert(ok_ == 1);
        Assert(0 <= cakePrice[ci]); Assert(cakePrice[ci] <= MaxPrice); }
    for (int ii = 0; ii < I; ++ii) { ok_ = scanf("%d", &ingPrice[ii]); Assert(ok_ == 1);
        Assert(0 <= ingPrice[ii]); Assert(ingPrice[ii] <= MaxPrice); }
    for (int ti = 0; ti < T; ++ti) { ok_ = scanf("%d", &toolPrice[ti]); Assert(ok_ == 1);
        Assert(0 <= toolPrice[ti]); Assert(toolPrice[ti] <= MaxPrice); }
    for (int ci = 0; ci < C; ++ci) for (int ii = 0; ii < I; ++ii) { 
        int amount; ok_ = scanf("%d", &amount); Assert(ok_ == 1);
        Assert(0 <= amount); Assert(amount <= MaxAmount);
        myint cost = myint(ingPrice[ii]) * amount;
        if (cakePrice[ci] < cost) cakePrice[ci] = 0;
        else cakePrice[ci] -= int(cost); }
    for (int ti = 0; ti < T; ++ti) nNeededBy[ti] = 0;
    for (int ci = 0; ci < C; ++ci) 
    {
        ok_ = scanf("%d", &nToolsNeeded[ci]); Assert(ok_ == 1); 
        Assert(0 <= nToolsNeeded[ci]); Assert(nToolsNeeded[ci] <= T);
        for (int i = 0; i < nToolsNeeded[ci]; ++i) {
            int ti; ok_ = scanf("%d", &ti); Assert(ok_ == 1);
            Assert(1 <= ti); Assert(ti <= T);
            --ti; needLists[ci][i] = ti; 
            neededByLists[ti][nNeededBy[ti]++] = ci; }
    }
    //
    const int source = 0, sink = 1 + C + T, firstCake = 1, firstTool = 1 + C;
    const int nVerts = sink + 1;
    enum { MaxVerts = MaxC + MaxT + 2 };
    // seen[u] = iterNo iff it was added to the queue in this iteration
    // prev[u] = the predecessor of u on the shortest path discovered by BFS; aug[u] = the maximum that this path can be augmented by
    int queue[MaxVerts], seen[MaxVerts], prev[MaxVerts]; 
    myint cakeFlow[MaxC], toolFlow[MaxT], cakeToolFlow[MaxC][MaxT], aug[MaxVerts];
    for (int ci = 0; ci < C; ++ci) cakeFlow[ci] = 0;
    for (int ti = 0; ti < T; ++ti) toolFlow[ti] = 0;
    for (int ci = 0; ci < C; ++ci) for (int ti = 0; ti < T; ++ti) cakeToolFlow[ci][ti] = 0;
    for (int u = 0; u < nVerts; ++u) seen[u] = 0;
    //
    myint totalFlow = 0;
    for (int iterNo = 1; ; ++iterNo)
    {
        // Edmonds-Karp: use breadth-first search to find an augmenting path.
        int head = 0, tail = 0;
        queue[tail++] = source; seen[source] = iterNo; prev[source] = -1; aug[source] = MaxPrice + 1;
        while (head < tail && seen[sink] != iterNo)
        {
            int u = queue[head++];
            if (u == 0) for (int ci = 0; ci < C; ++ci) {
                // source -> cake
                int v = firstCake + ci; if (seen[v] == iterNo) continue;
                myint spare = cakePrice[ci] - cakeFlow[ci]; if (spare <= 0) continue;
                aug[v] = min(aug[u], spare); seen[v] = iterNo; prev[v] = u; queue[tail++] = v; }
            else if (firstCake <= u && u < firstCake + C) 
                // cake -> tool
                for (int tii = 0, ci = u - firstCake; tii < nToolsNeeded[ci]; ++tii) {
                    int ti = needLists[ci][tii]; int v = firstTool + ti; if (seen[v] == iterNo) continue;
                    aug[v] = aug[u]; seen[v] = iterNo; prev[v] = u; queue[tail++] = v; }
            else if (firstTool <= u && u < firstTool + T) 
            {
                int ti = u - firstTool;
                // tool -> cake
                for (int cii = 0; cii < nNeededBy[ti]; ++cii) {
                    int ci = neededByLists[ti][cii]; int v = firstCake + ci; if (seen[v] == iterNo) continue;
                    myint spare = cakeToolFlow[ci][ti]; if (spare <= 0) continue;
                    aug[v] = min(aug[u], spare); seen[v] = iterNo; prev[v] = u; queue[tail++] = v; }
                // tool -> sink
                if (seen[sink] != iterNo)
                    if (myint spare = toolPrice[ti] - toolFlow[ti]; spare > 0) {
                        aug[sink] = min(aug[u], spare); seen[sink] = iterNo; prev[sink] = u; queue[tail++] = sink; }
            }
            else { Assert(u == sink); break; }
        }
        if (seen[sink] != iterNo) break; // no augmenting path exists
        // Augment this path.
        myint augBy = aug[sink]; Assert(augBy > 0);
        for (int v = sink; v != source; )
        {
            int u = prev[v]; 
            if (u == source) { 
                // source -> cake
                Assert(firstCake <= v && v < firstCake + C); int ci = v - firstCake;
                cakeFlow[ci] += augBy; Assert(0 <= cakeFlow[ci] && cakeFlow[ci] <= cakePrice[ci]); }
            else if (firstCake <= u && u < firstCake + C) {
                // cake -> tool
                Assert(firstTool <= v && v < firstTool + T); int ci = u - firstCake, ti = v - firstTool;
                cakeToolFlow[ci][ti] += augBy; Assert(0 <= cakeToolFlow[ci][ti]); }
            else if (firstTool <= u && u < firstTool + T) 
            {
                int ti = u - firstTool;
                // tool -> sink
                if (v == sink) { toolFlow[ti] += augBy; Assert(0 <= toolFlow[ti] && toolFlow[ti] <= toolPrice[ti]); }
                // tool -> cake
                else {
                    Assert(firstCake <= v && v < firstCake + C); int ci = v - firstCake;
                    cakeToolFlow[ci][ti] -= augBy; Assert(0 <= cakeToolFlow[ci][ti]); }
            }
            v = prev[v];
        }
        totalFlow += augBy;
    }
    myint totalFlowT = 0; for (int ti = 0; ti < T; ++ti) totalFlowT += toolFlow[ti];
    myint totalFlowC = 0; for (int ci = 0; ci < C; ++ci) totalFlowC += cakeFlow[ci];
    Assert(totalFlowT == totalFlow);
    Assert(totalFlowC == totalFlow);
    // total flow = (profit of cakes we haven't made) + (price of tools we have bought)
    //            = (profit of all cakes) - (profit of cakes we have made) + (price of tools we have bought)
    //            = (profit of all cakes) - [(profit of cakes we have made) - (price of tools we have bought)]
    //            = (profit of all cakes) - (net profit)
    myint totalCakePrice = 0; for (int ci = 0; ci < C; ++ci) totalCakePrice += cakePrice[ci];
    myint netProfit = totalCakePrice - totalFlow; Assert(netProfit >= 0);
    printf("%lld\n", (long long) netProfit);
    return 0;
}