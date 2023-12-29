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
#include <chrono>
using namespace std;

#define Assert assert
//#define Assert(x) 

enum { MaxN = 300'000, MaxTime = 1'000'000'000, Inf = MaxTime + 1 };

struct Lecture 
{
    int ti; // time when it was added
    int ai, bi; // start/end time
    bool isIn = false; //int biEff; // biEff = (isIn ? bi : Inf)
    // The following is the solution of running the greedy algorithm
    // from this lecture to the end of its block.
    int solCount = -1, solEndTime = -1; 
    bool operator < (const Lecture &other) const { return ai < other.ai; }
    bool operator < (int aiOther) const { return ai < aiOther; }
};

bool operator < (int ai, const Lecture &other) { return ai < other.ai; }

struct Operation { int time, lectureTi; bool add; };
struct Query { int nVisits = 0, lastVisit = -1; };

vector<short> firstByA; // firstByA[a] is the index of the first lecture in 'lectures' whose a_i is > a

struct Block
{
    int N;
    vector<Lecture> lectures; // All the N lectures in this block, ordered by a_i.
    vector<Operation> operations;

    void Init(const vector<Lecture> &lecturesSrc, int startIdx, int endIdx)
    {
        lectures.clear(); N = endIdx - startIdx; lectures.resize(N);
        for (int i = startIdx; i < endIdx; ++i) {
            auto &L = lectures[i - startIdx]; L = lecturesSrc[i];
            L.isIn = false; }
        sort(lectures.begin(), lectures.end());
    }

    void InitFirstByA(int srcTimeMax)
    {
        firstByA.resize(srcTimeMax);
        for (int a = 0, i = 0; a < srcTimeMax; ++a) {
            while (i < N && lectures[i].ai <= a) ++i;
            firstByA[a] = i; }
    }

    // Returns the index, in 'lectures', of the lecture which has the least value of biEff
    // amongst all those whose ai is > prevVisit.  Returns -1 if there is no such lecture.
    int RangeQuery(int prevVisit) const
    {
        int idx = (prevVisit < 0) ? 0 : (prevVisit >= firstByA.size()) ? N : firstByA[prevVisit];
        if (idx >= N) return -1; else return idx;
    }
};

int main()
{
    //freopen("attendance99.in", "rt", stdin);
    //freopen("attendance\\empty.in", "rt", stdin);
    auto tmStart = chrono::system_clock::now();
    int n; int ok_ = scanf("%d", &n); Assert(ok_ == 1);
    Assert(1 <= n); Assert(n <= MaxN);
    vector<Lecture> lectures; 
    vector<int> operations;
    vector<Query> queries(n); // queries[i] = calculates the solution we have to output after the first i operations
    for (int i = 0, addTime = 0; i < n; ++i)
    {
        int ai; ok_ = scanf("%d", &ai); Assert(ok_ == 1);
        if (ai < 0) { // removal
            ai = -ai;
            Assert(1 <= ai); Assert(ai <= lectures.size());
            auto &L = lectures[ai - 1]; Assert(L.isIn); L.isIn = false;
            operations.emplace_back(-ai); }
        else { // addition
            int bi; ok_ = scanf("%d", &bi); Assert(ok_ == 1);
            Assert(0 <= ai); Assert(ai <= bi); Assert(bi <= MaxTime);
            Lecture L; L.ti = addTime++; L.ai = ai; L.bi = bi; 
            L.isIn = true; lectures.push_back(L); 
            operations.emplace_back(addTime); }
    }
    // Compress the times.  If several a_i's fall between two b_i's, move them to the second b_i.
    struct Endpoint { int time, isEnd, i; };
    const int nLectures = lectures.size();
    vector<Endpoint> endpoints; endpoints.reserve(2 * nLectures);
    for (int i = 0; i < nLectures; ++i) { 
        const auto &L = lectures[i]; endpoints.push_back(Endpoint{L.ai, 0, i}); endpoints.push_back(Endpoint{L.bi, 1, i}); }
    sort(endpoints.begin(), endpoints.end(), [] (const auto &x, const auto &y) { 
        return x.time != y.time ? x.time < y.time : x.isEnd != y.isEnd ? x.isEnd < y.isEnd : x.i < y.i; });
    for (int i = endpoints.size() - 2; i >= 0; --i) 
        if (auto &E = endpoints[i]; E.time < endpoints[i + 1].time && ! E.isEnd) E.time = endpoints[i + 1].time;
    int prevSrcTime = -1, destTime = -1;
    for (auto &E : endpoints) {
        if (E.time > prevSrcTime) { prevSrcTime = E.time; ++destTime; }
        auto &L = lectures[E.i]; (E.isEnd ? L.bi : L.ai) = destTime; }
    const int maxTime = destTime + 1;
    // Sort the lectures by b_i.
    sort(lectures.begin(), lectures.end(), [] (const Lecture &x, const Lecture &y) { return x.bi < y.bi; });
    // Divide the lectures into blocks.
    int blockSize = 1; while (blockSize * blockSize < nLectures) ++blockSize;
    int nBlocks = (nLectures + blockSize - 1) / blockSize;
    vector<Block> blocks(nBlocks);
    for (int b = 0; b < nBlocks; ++b) blocks[b].Init(lectures, b * blockSize, min((b + 1) * blockSize, nLectures));
    lectures.clear(); // we shouldn't use these any more
    vector<pair<int, int>> lecturesByTi(nLectures, {-1, -1});
    for (int b = 0; b < nBlocks; ++b) { const auto &B = blocks[b]; 
        for (int i = 0; i < B.N; ++i) lecturesByTi[B.lectures[i].ti] = {b, i}; }
    // Divide the operations into blocks.        
    for (int i = 0; i < n; ++i)
    {
        const int op = operations[i];
        const int ti = (op < 0 ? -op : op) - 1;
        const int blockNo = lecturesByTi[ti].first, lectureNo = lecturesByTi[ti].second;
        Block &B = blocks[blockNo]; B.operations.push_back({i, ti, (op > 0)});
    }
    // Process the blocks.
    for (auto &B : blocks)
    {
        int nextQuery = 0; B.InitFirstByA(maxTime);
        for (int iOp = 0; iOp <= B.operations.size(); ++iOp)
        {
            // Perform the queries that precede this operation.
            int nextQuery2 = (iOp == B.operations.size()) ? n : B.operations[iOp].time;
            while (nextQuery < nextQuery2)
            {
                auto &Q = queries[nextQuery++];
                int nextIdx = B.RangeQuery(Q.lastVisit);
                if (nextIdx < 0) continue;
                auto &next = B.lectures[nextIdx];
                if (next.solCount <= 0) continue;
                Q.nVisits += next.solCount; Q.lastVisit = next.solEndTime;
            }
            if (iOp >= B.operations.size()) continue;
            // Perform this operation.
            auto &curOp = B.operations[iOp];
            const int lectureNo = lecturesByTi[curOp.lectureTi].second; 
            Lecture &L = B.lectures[lectureNo];
            if (curOp.add) { Assert(! L.isIn); L.isIn = true; }
            else { Assert(L.isIn); L.isIn = false; }
            // Recalculate the solutions for B's lectures.
            int bMin = Inf, bWhere = -1;
            for (int i = B.N - 1; i >= 0; --i)
            {
                // If we visit at time q, what is the latest time of the next visit?
                // That is the minimum b_j among lectures whose a_j is > q.
                auto &LL = B.lectures[i]; 
                if (LL.bi >= bMin || ! LL.isIn) { 
                    if (bWhere < 0) LL.solCount = 0, LL.solEndTime = -1;
                    else LL.solCount = B.lectures[bWhere].solCount, LL.solEndTime = B.lectures[bWhere].solEndTime; 
                    continue; }
                bMin = LL.bi; bWhere = i;
                int nextIdx = B.RangeQuery(LL.bi);
                if (nextIdx < 0) { LL.solCount = 1; LL.solEndTime = LL.bi; }
                else {
                    Assert(nextIdx > i);
                    const auto &next = B.lectures[nextIdx]; 
                    LL.solCount = next.solCount + 1; LL.solEndTime = (next.solCount > 0) ? next.solEndTime : LL.bi; }
            }
        }
    }
    for (const auto &Q : queries) printf("%d\n", Q.nVisits);
    fprintf(stderr, "%.6f sec\n", chrono::duration<double, chrono::seconds::period>(chrono::system_clock::now() - tmStart).count());
    return 0;
}