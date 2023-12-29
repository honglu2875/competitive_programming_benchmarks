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

enum { MaxN = 30'000, MaxQ = 300'000, MaxL = 300'000, MaxT = 1'000'000'000 };

int main()
{
    //freopen("laundry99.in", "rt", stdin);
    //freopen("laundry\\random_small4.in", "rt", stdin);
    typedef long long int llint;
    typedef uint_fast64_t ullint; enum { ULL_BITS = 64 };
    struct Sheet { llint d, a, b, c, D; }; // a_i = t^fast_i, b_i = t^slow_i
        // Sheet.c is the maximum of this.b and of the values of a for all sheets whose b is > this.b.
        // Sheet.D is the sum of d for all sheets whose b is > this.b.
    struct Query { llint L, r = -1; }; // L_j and the corresponding result
    // Read the input data.
    llint n, q, dSum = 0, ok_ = scanf("%lld %lld", &n, &q); Assert(ok_ == 2);
    Assert(1 <= n); Assert(n <= MaxN); Assert(1 <= q); Assert(q <= MaxQ);
    vector<Sheet> sheets(n);
    for (auto &S : sheets) { 
        ok_ = scanf("%lld %lld %lld", &S.d, &S.a, &S.b); Assert(ok_ == 3);
        Assert(1 <= S.d); Assert(S.d <= MaxL);
        Assert(1 <= S.a); Assert(S.a <= S.b); Assert(S.b <= MaxT); 
        dSum += S.d; }
    vector<Query> queries(q); llint lMax = 0;
    for (auto &Q : queries) { 
        ok_ = scanf("%lld", &Q.L); Assert(ok_ == 1);
        Assert(1 <= Q.L); Assert(Q.L <= MaxL);
        lMax = max(lMax, Q.L); }
    // Sort the sheets by b_i.
    sort(sheets.begin(), sheets.end(), [] (const auto &x, const auto &y) { return x.b < y.b; });
    // Suppose we want to dry everything in at most T units of time.
    // Any sheet that has b_i > T will have to be hung over two lines;
    // but sheets with b_i <= T can be hung over one line and there can't
    // be any advantage in hanging them over two lines.  Hence the optimal
    // solution is always of the following form: the first k sheets (for 
    // some k from 1..n; remember that the sheets are now ordered by b_i)
    // will be hung over one line, the others in two.  In this case the
    // time needed for all the sheets to dry will be
    //   c_k = max { b_1, ..., b_k, a_{k+1}, ..., a_n }
    //       = max { b_k, a_{k+1}, ..., a_n }.
    // We'll store this in sheets[k].c.  Note that c_k <= c_{k+1}: 
    // if c_k = b_k, then c_{k+1} >= b_{k+1} >= b_k = c_k;
    // if c_k = a_{k+1}, then c_{k+1} >= b_{k+1} >= a_{k+1} = c_k;
    // and if c_k = a_i for some i > k + 1, then c_{k+1} >= a_i = c_k.
    // Thus, the candidate solutions c_k form an increasing sequence.
    llint aMax = 0;
    for (llint i = n - 1, D = 0; i >= 0; --i) { auto &S = sheets[i];
        S.c = aMax = max(aMax, S.a);
        S.D = D; D += S.d; }
    for (int i = 0; i < n; ++i) { auto &S = sheets[i];
        S.c = max(S.c, S.b);
        if (i > 0) Assert(S.c >= sheets[i - 1].c); }
    // But when is c_k a viable candidate for a given line length L?
    // Let E_k := d_1 + ... + d_k and D_k := d_{k+1} + ... + d_n.
    // [Note that D_k is already ready in sheets[k].D.]
    // Thus the sheets k+1..n, which use both lines, cover D_k of each line
    // and the sheets 1..k have L - D_k space on each line.
    // Let f(d, k) be a boolean value which indicates whether some 
    // subset of sheets 1..k exists such that their d_i sum up to exactly d.
    // [This is a backpack problem that can be solved with dynamic programming:
    // f(d, k) = f(d, k - 1) or f(d - d_i, k - 1).]
    // When asking whether the sheets 1..k can be dried on two lines of length L - D_k,
    // we can find the maximum u_k <= L - D_k for which f(u_k, k) is true;
    // these sheets can be placed on the first line, and we then just have to check
    // if the remaining sheets have enough space on the second line, i.e. if E_k - u_k <= L - D_k.
    // We can rewrite this condition as E_k + D_k <= L + u_k.  Note that the left side
    // is now constant, simply the sum of the d_i of all sheets.  What happens to the 
    // right side if we slowly increase L?  Remember that u_k = max { d <= L - D_k : f(d, k) }.
    // So when L increases, additional d's come into consideration for this max {...}, and
    // consequenly u_k can only increase, never decrease.  Thus as L increaes, the right side
    // of our inequality E_k + D_k <= L + u_k also increases.  This means that if
    // candidate c_k is viable for a given L, it is also viable for every greater L.
    // Let L_k be the minimum L for which candidate c_k is viable.  
    int fSize = lMax / ULL_BITS + 1; // we want bits f[0] ... f[lMax] to be available
    vector<ullint> f(fSize, 0); f[0] = 1; llint dMax = 0;
    //vector<bool> fTest_(lMax + 1, false); fTest_[0] = true;
    vector<int> highestBit(1 << 16, -1); for (int k = 0; k < 16; ++k) for (int i = 0; i < (1 << k); ++i) highestBit[i | (1 << k)] = k;
    struct Cand { llint c, L; }; // c_k and L_k for some value of k
    vector<Cand> cands;
    // First we have a candidate where all the sheets use both lines.
    if (dSum <= lMax) cands.push_back({aMax, dSum});
    // Next consider candidates where the some sheets use one line.
    for (int k = 0; k < n; ++k)
    {
        auto &S = sheets[k];
        // dMax is the value above which all f(d, k) are guaranteed to be false.
        dMax = min(lMax, dMax + S.d);
        // Currently f[d] = f(d, k-1).  Calculate f(d, k) and store them into f[d].
        // We want to perform the operation f[d] |= f[d - d_k] for each bit d.
        // But our f is actually an array of ullint's.  
        int ofs = int(S.d + ULL_BITS - 1) / ULL_BITS;
        int shift1 = (ULL_BITS - (S.d % ULL_BITS)) % ULL_BITS;
        int shift2 = ULL_BITS - shift1;
        for (int i = fSize - 1; i >= ofs - 1; --i)
        {
            //if (k == 9 && i == 31) printf("!");
            f[i] |= (i >= ofs ? f[i - ofs] >> shift1 : 0) | (shift1 ? f[i - ofs + 1] << shift2 : 0); // doing just 'x << shift2' when shift2 == 64 might give you the result x instead of the expected 0
            //if (k == 9 && i == 31) printf("!");
        }
        /*
        for (llint d = dMax; d >= S.d; --d) {
            if (fTest_[d - S.d]) fTest_[d] = true;
            Assert((fTest_[d] ? 1 : 0) == ((f[d / ULL_BITS] >> (d % ULL_BITS)) & 1)); } 
        */
        // Find the least L such that E_k + D_k <= L + u_k, where u_k = max { d <= L - D_k : f(d, k) }.
        // (Obviously we should only consider L >= D_k, otherwise the lines are
        // too short even just for the sheets that require both lines.)
        // This L will be the L_k we're looking for.
        // Once we have c_k and L_k for every candidate candidate: when a query L comes,
        // the solution is min_k { c_k : c_k is a viable candidate }
        //               = min_k { c_k : L >= L_k },
        // but since c_1 <= c_2 <= ... <= c_n, the minimum will always be achieved by 
        // the first viable candidate, i.e. the one with the least k.  So the answer is c_k
        // where k = min { k : L >= L_k }.  If we have two candidates k and k' where
        // k < k' and L_k <= L_k', then any L which is >= L_k' is also >= L_k and
        // will always prefer candidate k over k'; hence k' may be discarded.
        // In our case this means (since we're generating candidates in increasing order
        // of k) that the last candidate currently in 'cands' must be greater than the
        // current one, otherwise there's no point in keeping the current one.
        // This means that our loop over L can stop when it reaches the L_k of the previous candidate.
        llint Lk = -1;
        llint lFrom = S.D, lTo = cands.empty() ? lMax : cands.back().L;
        if (lTo < lFrom) continue;
        // So in principle our loop can be:
        //     for (L = lFrom, u_k = 0; L <= lTo; ++L) {
        //       if f(L - D_k, k) u_k = L - D_k;
        //       if (E_k + D_k <= L + u_k) { L_k = L; break; } }
        // But since our function f(_, k) is packed into a vector of 64-bit integers,
        // it's easier to use L - D_k as a counter:
        //     for (LL = 0, u_k = 0; LL <= lTo - S.D; ++L) {
        //       if (f(LL, k) u_k = LL;
        //       if (E_k + D_k <= LL + D_k + u_k) { L_k = LL + D_k; break; } }
        // In the last line, we can simplify the condition to 'E_k <= LL + u_k'.
        llint llFrom = lFrom - S.D, llTo = lTo - S.D; 
        llint Ek = dSum - S.D; 
        llint fFrom = llFrom / ULL_BITS, fTo = (llTo + 1) / ULL_BITS;
        Assert(llFrom == 0); Assert(fFrom == 0);
        for (llint fi = fFrom, uk = 0; fi <= fTo; ++fi) 
        {
            ullint fCur = f[fi];
            llint ukNext;
            if (fCur & 0xffff'0000'0000'0000ULL) ukNext = fi * ULL_BITS + 48 + highestBit[(fCur >> 48) & 0xffff];
            else if (fCur & 0xffff'0000'0000ULL) ukNext = fi * ULL_BITS + 32 + highestBit[(fCur >> 32) & 0xffff];
            else if (fCur & 0xffff'0000ULL) ukNext = fi * ULL_BITS + 16 + highestBit[(fCur >> 16) & 0xffff];
            else if (fCur) ukNext = fi * ULL_BITS + highestBit[fCur]; 
            else ukNext = uk;
            // The current entry of 'f', namely f[fi], represents the range of LL's from fi * ULL_BITS to (fi + 1) * ULL_BITS + 1.
            // During this time, u_k will grow from its current value 'uk' to 'ukNext'.
            // Perhaps even that won't be enough to satisfy the condition E_k <= LL + u_k.
            if (Ek > (fi + 1) * ULL_BITS + ukNext) { uk = ukNext; continue; }
            // Otherwise go bit by bit.
            for (llint LL = fi * ULL_BITS, bit = 0; bit < ULL_BITS; ++bit, ++LL, fCur >>= 1)  {
                if (fCur & 1) uk = LL;
                if constexpr (false) if (k == 16 && LL + S.D >= 299930) printf("LL = %lld, L = %lld, uk = %lld;  Ek = %lld, LL + uk = %lld\n", LL, LL + S.D, uk, Ek, LL + uk);
                if (Ek <= LL + uk) { Lk = LL + S.D; break; }
            }
            if (Lk >= 0) break;
            Assert(uk == ukNext);
        }
        if constexpr (false) 
        {
            llint Lk_ = -1;
            for (llint L = lFrom, uk = 0; L <= lTo; ++L) {
                //if (fTest_[L - S.D]) uk = L - S.D;
                if ((f[(L - S.D) / ULL_BITS] >> ((L - S.D) % ULL_BITS)) & 1) uk = L - S.D;
                if constexpr (false) if (k == 16 && L >= 299930) printf("or LL = %lld, L = %lld, uk = %lld; dSum = %lld, L + uk = %lld\n", L - S.D, L, uk, dSum, L + uk);
                if (dSum <= L + uk) { Lk_ = L; break; } }
            Assert(Lk_ == Lk);
            //fprintf(stderr, "k = %d, Lk = %lld, Lk_ = %lld\n", k, Lk, Lk_);
        }
        if constexpr (false) if (k % 100 == 0) fprintf(stderr, "k = %d -> Lk = %lld    \r", k, Lk);
        //fprintf(stderr, "k = %d -> ck = %lld, Lk = %lld\n", k, S.c, Lk);
        if (Lk < 0) continue;
        //if (! cands.empty() && cands.back().L <= Lk) continue;  // this does happen sometimes (tested)
        // This will leave us with a list of candidates where as k increases, c_k increases and L_k decreases.
        cands.push_back({S.c, Lk});
    }
    // For each query, we can find the best candidate using binary search.
    for (auto &Q : queries)
    {
        // We want the first candidate whose L_k is <= Q.L.
        auto it = lower_bound(cands.begin(), cands.end(), Q.L,
            [] (const Cand &x, llint y) { return x.L > y; });
        printf("%lld\n", it == cands.end() ? -1 : it->c);
    }
    return 0;
}
