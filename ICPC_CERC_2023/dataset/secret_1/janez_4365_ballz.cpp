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
#include <map>
using namespace std;

#define Assert assert
//#define Assert(x) 

enum { MaxN = 50, MaxCoord = 10000 };
struct Point { int x, y; };

int CCW(Point A, Point B, Point C)
{
    return (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
}

long double SqrtDist(Point A, Point B)
{
    int dx = A.x - B.x, dy = A.y - B.y; return sqrt((long double) (dx * dx + dy * dy));
}

int main()
{
    // Read the input data.
    int n, ok_; 
    string line; getline(cin, line);
    ok_ = sscanf(line.c_str(), "%d", &n); Assert(ok_ == 1); Assert(2 <= n); Assert(n <= MaxN);
    string sex; getline(cin, sex); Assert(sex.length() == n);
    for (char c : sex) Assert(c == 'B' || c == 'G');
    vector<Point> Ts(n), bySex[2]; 
    for (int i = 0; i < n; ++i){
        line.clear(); getline(cin, line);
        auto &T = Ts[i];
        ok_ = sscanf(line.c_str(), "%d %d", &T.x, &T.y); Assert(ok_ == 2);
        Assert(-MaxCoord <= T.x); Assert(T.x <= MaxCoord);
        Assert(-MaxCoord <= T.y); Assert(T.y <= MaxCoord);
        bySex[sex[i] == 'B' ? 0 : 1].push_back(T);
    }
    line.clear(); getline(cin, line); Assert(line.empty()); Assert(cin.eof());
    // Check for strict convexity.
    for (int i = 0; i < n; ++i)
    {
        int side = 0;
        for (int d = 2; d < n; ++d)
        {
            int ccw = CCW(Ts[i], Ts[(i + 1) % n], Ts[(i + d) % n]);
            Assert(ccw != 0);
            if (ccw > 0) ccw = 1; else ccw = -1;
            if (d == 2) side = ccw; else Assert(ccw == side);
        }
    }
    // Calculate and output the result.
    long double result = 0;
    for (int Sex = 0; Sex < 2; ++Sex)
    {
        auto &v = bySex[Sex]; int m = v.size(); 
        Assert(m % 2 == 0); m /= 2;
        for (int i = 0; i < m; ++i)
            result += SqrtDist(v[i], v[i + m]);
    }
    /*
    // Ballz deep!
    constexpr long double RelErr = 1.000'001, AbsErr = 0.000'001;
    long double minResult = min(result - AbsErr, min(result * RelErr, result / RelErr));
    long double maxResult = max(result + AbsErr, max(result * RelErr, result / RelErr));
    mt19937_64 r(69'105'420);
    long double c = uniform_real_distribution<long double>(0, 1)(r);
    c = 1; 
    c = 0.001 + 0.998 * c;
    result = minResult * c + maxResult * (1 - c);
    */
    printf("%.12Lf\n", result); return 0;
}