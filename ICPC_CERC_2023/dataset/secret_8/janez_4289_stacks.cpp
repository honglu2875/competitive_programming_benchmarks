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
#include <unordered_set>
using namespace std;

#define Assert assert
//#define Assert(x) 

struct Case
{
    int n, k;
    vector<int> C; // capacities
    vector<vector<int>> S; // stacks; only nonzero values are stored
    struct Where { int stack, pos; };
    vector<Where> where;
    // 'nonfullStacks' is a list of non-full stacks, in arbitrary order.  whereInNonfullStacks[s] is -1 if stack 's' is full,
    // otherwise it tells you that nonfullStacks[whereInNonfullStacks[s]] == s.
    vector<int> nonfullStacks, whereInNonfullStacks;

    void Init(int n_, int k_, vector<int> capacities, vector<vector<int>> cardsZeroBased)
    {
        n = n_; k = k_; 
        C = move(capacities); S = move(cardsZeroBased);
        Assert(C.size() == k); Assert(S.size() == k);
        int cSum = 0, cMax = 0; 
        where.clear(); where.resize(n, {-1, -1});
        int nSeen = 0;
        nonfullStacks.clear(); whereInNonfullStacks.clear();
        whereInNonfullStacks.resize(k, -1); 
        for (int i = 0; i < k; ++i) 
        {
            Assert(C[i] > 0); Assert(S[i].size() <= C[i]);
            cSum += C[i]; cMax = max(cMax, C[i]);
            for (int j = 0; j < S[i].size(); ++j) {
                int x = S[i][j]; auto &Wx = where[x];
                Assert(Wx.stack < 0); Wx.stack = i; Wx.pos = j; ++nSeen; }
            if (S[i].size() < C[i]) { whereInNonfullStacks[i] = nonfullStacks.size(); nonfullStacks.emplace_back(i); }
        }
        Assert(nSeen == n);
        Assert(n + cMax <= cSum);
    }

    void Read(FILE *f = nullptr)
    {
        if (! f) f = stdin;
        int ok_ = fscanf(f, "%d %d", &n, &k); Assert(ok_ == 2);
        C.resize(k); S.resize(k);
        int cSum = 0, cMax = 0;
        where.clear(); where.resize(n, {-1, -1});
        int nSeen = 0;
        nonfullStacks.clear(); whereInNonfullStacks.clear();
        whereInNonfullStacks.resize(k, -1); 
        for (int i = 0; i < k; ++i)
        {
            ok_ = fscanf(f, "%d", &C[i]); Assert(ok_ == 1);
            Assert(C[i] > 0); S[i].resize(C[i]);
            cSum += C[i]; cMax = max(cMax, C[i]);
            bool ended = false;
            for (int j = 0; j < C[i]; ++j) 
            {
                int &Sij = S[i][j];
                ok_ = fscanf(f, "%d", &Sij); Assert(ok_ == 1);
                Assert(Sij >= 0); Assert(Sij <= n); --Sij;
                if (Sij < 0) ended = true;
                else 
                { 
                    Assert(! ended); // zeros must be at the top
                    auto &W = where[Sij]; Assert(W.stack < 0); // each card must appear exactly once
                    W.stack = i; W.pos = j; ++nSeen; 
                } 
            }
            while (! S[i].empty() && S[i].back() < 0) S[i].pop_back();
            if (S[i].size() < C[i]) { whereInNonfullStacks[i] = nonfullStacks.size(); nonfullStacks.emplace_back(i); }
        }
        Assert(nSeen == n);
        Assert(n + cMax <= cSum);
        int dummy; ok_ = fscanf(f, "%d", &dummy); Assert(ok_ <= 0);
    }

    int nMoves = 0; FILE *fOut = nullptr;
    void Move(int from, int to) 
    {
        Assert(0 <= from); Assert(from < k); Assert(0 <= to); Assert(to < k);
        auto &S1 = S[from], &S2 = S[to];
        Assert(! S1.empty());
        if (from == to) return;
        Assert(S2.size() < C[to]);
        if (S1.size() >= C[from]) {
            // S1 was full but won't be any more, so add it to the nonfull stacks list.
            Assert(whereInNonfullStacks[from] == -1);
            whereInNonfullStacks[from] = nonfullStacks.size(); nonfullStacks.emplace_back(from); }
        int x = S1.back(); S1.pop_back();
        Assert(where[x].stack == from); Assert(where[x].pos == S1.size());
        where[x].stack = to; where[x].pos = S2.size();
        S2.emplace_back(x); ++nMoves;
        if (S2.size() >= C[to]) { 
            // S2 has become full with this move, so remove it from the nonfull stacks list.
            int ei = whereInNonfullStacks[to]; Assert(ei >= 0); Assert(nonfullStacks[ei] == to); 
            int other = nonfullStacks.back(); nonfullStacks[ei] = other; whereInNonfullStacks[other] = ei;
            nonfullStacks.pop_back(); whereInNonfullStacks[to] = -1; }
        if (fOut) fprintf(fOut, "%d %d\n", from + 1, to + 1);
    }

    // Finds a stack with some empty space (i.e. a non-full stack).
    int FindEmpty(int avoid1 = -1, int avoid2 = -1) const
    {
        for (int stack : nonfullStacks) if (stack != avoid1 && stack != avoid2) return stack;
        return -1;
    }

    bool IsOnTop(int x) const { return where[x].pos == S[where[x].stack].size() - 1; }
    bool IsOnTop(int x, int stack) const { return where[x].stack == stack && where[x].pos == S[stack].size() - 1; }
    bool IsFull(int stack) const { return S[stack].size() >= C[stack]; }

    int Solve(FILE *fOut_)
    {
        fOut = fOut_; nMoves = 0;
        for (int destStack = 0, x = 0; destStack < k && x < n; ++destStack) for (int destPos = 0; destPos < C[destStack] && x < n; ++destPos, ++x)
        {
            // Earlier stacks, and earlier spots on the current stack, already contain 
            // cards from 0 to x - 1 in their desired positions, so 'x' certainly can't be there.
            auto &Wx = where[x];
            Assert(Wx.stack > destStack || (Wx.stack == destStack && Wx.pos >= destPos));
            // Perhaps 'x' already is where it should be.
            if (Wx.stack == destStack && Wx.pos == destPos) continue;
            // Ensure that 'x' is at the top of its stack, by moving cards above it
            // to any other stacks where space is available.
            while (! IsOnTop(x)) Move(Wx.stack, FindEmpty(Wx.stack));
            // Move 'x' to the top of the correct stack.
            if (Wx.stack != destStack) 
            {
                // If the destination stack is full, move one card from it to any other stack (except
                // the one that 'x' is currently on top of, of course).
                if (IsFull(destStack)) Move(destStack, FindEmpty(Wx.stack));
                Assert(IsOnTop(x));
                // Move x to the destination stack.
                Move(Wx.stack, destStack);
            }
            // Now 'x' is on top of the correct stack, but there may still be unwanted cards below it.
            Assert(IsOnTop(x, destStack));
            while (Wx.pos > destPos)
            {
                Assert(S[destStack].size() >= 2);
                // 'x' is on the correct stack but there is at least one unwanted card on it.
                // This means that 'destStack' contains at least two cards, so there must be at least
                // two empty spots in the other stacks.  
                int s1 = FindEmpty(destStack); Assert(s1 > destStack);
                int s2 = FindEmpty(destStack, s1);
                // It could happen, however, that all the empty spots are concentrated in a single stack.
                if (s2 < 0)
                {
                    // All the empty spots are in 's1' (and perhaps on 'destStack'; but that one also has
                    // at least two cards), which therefore has at least two empty spots.
                    if (nonfullStacks.size() == 1) Assert(nonfullStacks[0] == s1);
                    else { Assert((nonfullStacks[0] == s1 && nonfullStacks[1] == destStack) || (nonfullStacks[1] == s1 && nonfullStacks[0] == destStack)); }
                    Assert(S[s1].size() <= C[s1] + 2);
                    // Find some full stack and move one card from it to 's1'.
                    s2 = destStack + 1; if (s2 == s1) ++s2;
                    if (s2 < k) 
                    {
                        Assert(IsFull(s2)); 
                        Move(s2, s1); 
                    }
                    else {
                        // But wait!  It could be the case that 'destStack' is the last stack but one (k - 2),
                        // 's1' is the last stack (k - 1) and all earlier stacks are full.  In this case we
                        // can use one of the earlier stacks temporarily as 's2', but must be careful to
                        // restore it to its previous condition later.
                        s2 = 0; Assert(IsFull(s2)); Assert(s2 < destStack);
                        Move(s2, s1);
                    }
                }
                Assert(s1 > destStack); Assert(s2 != destStack); Assert(s1 != s2);
                Assert(! IsFull(s1)); Assert(! IsFull(s2));
                // Move 'x' temporarily to 's2', then move the card that used to be just below it to 's1',
                // then move 'x' back to the correct stack.
                Move(destStack, s2); Assert(IsOnTop(x, s2));
                Move(destStack, s1); 
                Move(s2, destStack); Assert(IsOnTop(x, destStack));
                // Now 's2' is in the same condition as it was before we moved 'x' temporarily onto it.
                // But perhaps 's2' is one of those earlier stacks that were already in their final state
                // before we moved one card from it to 's1', in which case we must now undo this.
                if (s2 < destStack)
                {
                    Assert(! IsFull(s2));
                    // But wait!  The card which we moved from 's2' to 's1' is now buried under the
                    // card that we just moved from 'destStack' (the one that used to be just below 'x').
                    // So we have to temporarily move that card back to 'destStack' so as to expose the
                    // one that needs to be moved back to 's2'.
                    Move(s1, destStack);
                    Move(s1, s2);
                    Move(destStack, s1);
                }
            }
        }
        for (int destStack = 0, x = 0; destStack < k && x < n; ++destStack) for (int destPos = 0; destPos < C[destStack]; ++destPos, ++x)
            if (x < n) Assert(S[destStack][destPos] == x);
            else Assert(S[destStack].size() <= destPos);
        if (fOut) fprintf(fOut, "0 0\n");
        return nMoves;
    }
};

void GenRndCards(mt19937_64 &r, int n, const vector<int> &C, vector<vector<int>> &stacks)
{
    int cSum = 0, cMax = 0; for (int c : C) cSum += c, cMax = max(cMax, c); 
    Assert(n <= cSum - cMax);
    vector<int> cards(cSum); for (int i = 0; i < cSum; ++i) cards[i] = i;
    shuffle(cards.begin(), cards.end(), r);
    const int k = C.size();
    stacks.resize(k); for (auto &s : stacks) s.clear();
    for (int i = 0, x = 0; i < k; ++i) for (int j = 0; j < C[i]; ++j, ++x)
        if (cards[x] < n) stacks[i].emplace_back(cards[x]); 
}

int TestRandom(int argc, char **argv)
{
    vector<int> maxMoves(100, 0), maxWhen(100, -1);
    for (int nIter = 0; ; ++nIter)
    {
        mt19937_64 r(567890 + nIter);
        int n = uniform_int_distribution(1, 20)(r);
        int k = uniform_int_distribution(3, 10)(r);
        n = 18; k = 3; // it seems that the largest number of moves is required when k = 3
        if (argc > 1) sscanf(argv[1], "%d", &n); 
        vector<int> C(k, 1); int cSum = k, cMax = 1;
        while (n > cSum - cMax) {
            int i = uniform_int_distribution(0, k - 1)(r); C[i] += 1; 
            ++cSum; cMax = max(cMax, C[i]); }
        int nExtra = uniform_int_distribution(0, 5)(r);
        nExtra = 0; // for the maximum number of moves
        for (int u = 0; u < nExtra; ++u) {
            int i = uniform_int_distribution(0, k - 1)(r); C[i] += 1; ++cSum; }
        //
        vector<vector<int>> stacks; GenRndCards(r, n, C, stacks);
        //
        Case c; c.Init(n, k, C, stacks);
        int nMoves = c.Solve(nullptr);
        if (nMoves > maxMoves[n]) { maxMoves[n] = nMoves; maxWhen[n] = nIter; 
            printf("%d]  n = %d, k = %d, cSum = %d (%d extra) [%d, %d, %d] -> %d moves\n", nIter, n, k, cSum, nExtra, C[0], C[1], C[2], nMoves); }
    }
    return 0;
    /*
    245606]  n = 15, k = 3, cSum = 29 (0 extra) [2, 14, 13] -> 478 moves
    624419]  n = 16, k = 3, cSum = 30 (0 extra) [2, 14, 14] -> 548 moves
    587941]  n = 17, k = 3, cSum = 31 (0 extra) [3, 14, 14] -> 586 moves
    213748]  n = 18, k = 3, cSum = 33 (0 extra) [4, 15, 14] -> 656 moves
    109352]  n = 19, k = 3, cSum = 38 (0 extra) [3, 19, 16] -> 712 moves
    11802]   n = 20, k = 3, cSum = 38 (0 extra) [3, 18, 17] -> 838 moves
    1343232] n = 25, k = 3, cSum = 46 (0 extra) [5, 20, 21] -> 1163 moves
    1548621]  n = 30, k = 3, cSum = 54 (0 extra) [6, 24, 24] -> 1648 moves

    TestRandom2:
    iter 29664] n = 20, k = 3, cSum = 40, cMax = 20 [1, 20, 20] -> 1136 moves
    iter 12919] n = 25, k = 3, cSum = 50, cMax = 25 [1, 25, 25] -> 1758 moves
    iter 6162] n = 30, k = 3, cSum = 60, cMax = 30 [1, 30, 30] -> 2477 moves

    iter 58769] n = 20, k = 3, cSum = 41, cMax = 20 [1, 20, 20] -> 1204 moves  // random
    iter 0] n = 20, k = 3, cSum = 41, cMax = 20 [1, 20, 20] -> 1247 moves      // initial stacks: [], [], [0, ..., n-1]  (i.e. 1 at the bottom)
    iter 2751409] n = 20, k = 3, cSum = 41, cMax = 20 [1, 20, 20] -> 1246 moves //  [ ] [ 16 13 19 14 0 17 12 15 11 18 ] [ 2 1 6 3 5 4 8 7 9 10 ]
    iter 76103] n = 25, k = 3, cSum = 51, cMax = 25 [1, 25, 25] -> 1806 moves
    
    iter 70544] n = 30, k = 3, cSum = 61, cMax = 30 [1, 30, 30] -> 2521 moves // [ 28 ] [ 23 16 27 25 19 26 29 24 0 17 18 ] [ 4 5 6 3 10 15 21 7 1 12 11 13 2 20 9 22 14 8 ]
    iter 332829] n = 30, k = 3, cSum = 61, cMax = 30 [1, 30, 30] -> 2544 moves //  [ ] [ 21 15 20 0 29 28 16 27 11 24 19 25 23 17 10 26 ] [ 5 1 18 2 22 3 12 6 7 9 4 8 13 14 ]
    iter 15999280] n = 30, k = 3, cSum = 61, cMax = 30 [1, 30, 30] -> 2640 moves // [ 11 ] [ 29 22 25 24 28 23 7 8 27 26 20 18 21 17 ] [ 2 4 19 1 5 6 3 0 10 15 9 16 14 12 13 ]
    iter 17792652] n = 30, k = 3, cSum = 61, cMax = 30 [1, 30, 30] -> 2682 moves // [ 18 ] [ 29 25 23 22 28 24 26 17 27 20 21 19 ] [ 10 2 12 14 1 6 7 4 8 15 11 5 16 13 3 0 9 ]

    iter 0] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 10 moves // [], [], [9..0]
    iter 0] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 84 moves  // [], [0..9], []
    iter 0] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 237 moves   // [], [], [0..9]
    iter 0] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 289 moves  // [], [9..0], []
    iter 15710] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 300 moves  // [ 8 ] [ 9 0 7 6 ] [ 1 2 3 4 5 ]
    iter 619226] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 301 moves  // [ ] [ 9 8 0 6 5 4 7 ] [ 1 2 3 ]
    iter 1579461] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 303 moves // [ 9 ] [ 0 8 7 6 5 ] [ 1 2 3 4 ]
    iter 40005679] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 305 moves  // [ ] [ 0 8 7 6 5 4 9 ] [ 1 2 3 ]
    iter 126786890] n = 10, k = 3, cSum = 21, cMax = 10 [1, 10, 10] -> 307 moves //  [ ] [ 0 8 7 6 5 4 3 2 9 ] [ 1 ]
    */
}

int TestRandom2(int argc, char **argv)
{
    int n = 20; if (argc > 1) sscanf(argv[1], "%d", &n); 
    int maxMoves = 0;
    for (int nIter = 0; ; ++nIter)
    {
        mt19937_64 r(654 + nIter);
        vector<int> C(3);
        /*
        for (C[0] = 1; C[0] <= 3; ++C[0])
        for (C[1] = 1; C[1] <= C[0] + n; ++C[1])
        {
            int cSum = C[0] + C[1], cMax = max(C[0], C[1]);
            C[2] = 1;
            while (n > cSum - cMax && C[2] <= max(C[0], C[1])) { 
                ++C[2]; cMax = max(cMax, C[2]); ++cSum; }
            C[0] = (n + 1) / 2; C[1] =C[0]; C[2] = C[0]; cSum = 3 * C[0]; cMax = C[0]; Assert(n <= cSum - cMax);
            if (n > cSum - cMax) continue;
            vector<vector<int>> stacks; GenRndCards(r, n, C, stacks);
            Case c; c.Init(n, 3, C, stacks);
            int nMoves = c.Solve(nullptr);
            if (nMoves > maxMoves) { 
                maxMoves = nMoves;
                printf("iter %d] n = %d, k = %d, cSum = %d, cMax = %d [%d, %d, %d] -> %d moves\n",
                    nIter, n, 3, cSum, cMax, C[0], C[1], C[2], nMoves);
            }
        }
        */
        // Experiments show that for a given n, the largest number of moves is required when
        // there are just k = 3 stacks; and where their capacities are [1, n, n].  Allowing higher
        // capacities, or having the capacity 1 somewhere else than on the first stack, 
        // reduces the number of moves required.
        C[0] = n; C[1] = n; C[2] = n;
        //C[uniform_int_distribution(0, 2)(r)] = 1;
        C[0] = 1;
        int cSum = C[0] + C[1] + C[2], cMax = max(C[0], max(C[1], C[2]));
        Assert(n <= cSum - cMax);
        vector<vector<int>> stacks; GenRndCards(r, n, C, stacks);
        //stacks.clear(); stacks.resize(3); for (int x = 0; x < n; ++x) stacks[2].emplace_back(x);
        //stacks.clear(); stacks.resize(3); for (int x = 0; x < n; ++x) stacks[2].emplace_back(n - 1 - x);
        Case c; c.Init(n, 3, C, stacks);
        int nMoves = c.Solve(nullptr);
        if (nMoves > maxMoves) { 
            maxMoves = nMoves;
            printf("iter %d] n = %d, k = %d, cSum = %d, cMax = %d [%d, %d, %d] -> %d moves\n",
                nIter, n, 3, cSum, cMax, C[0], C[1], C[2], nMoves);
            for (auto &S : stacks) { printf(" ["); for (int x : S) printf(" %d", x); printf(" ]"); }
            printf("\n");
        }
    }
}

int main(int argc, char** argv)
{
    if (false) return TestRandom(argc, argv);
    if (false) return TestRandom2(argc, argv);
    //freopen("stacks99.in", "rt", stdin);
    Case c; c.Read();
    int nMoves = c.Solve(stdout);
    fprintf(stderr, "n = %d, k = %d; solved in %d moves\n", c.n, c.k, nMoves);
    return 0;
}