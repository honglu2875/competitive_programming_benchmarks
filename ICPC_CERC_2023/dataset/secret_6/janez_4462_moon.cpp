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

//#define Assert assert
#define Assert(x) 

enum { MaxT = 1000, MaxCoord = 1000, MaxR = 1000 };
enum { verbose = false };
//constexpr double pi = 3.14159265358979;
const long double pi = 4 * atan2((long double) 1, (long double) 1);
enum { MUST_TOUCH_EDGE = false };

namespace Tomaz
{

    typedef long long int64;
    typedef pair<int,int> PII;
    typedef vector<int> VI;
    typedef array<int,3> III;
    typedef vector<string> VS;

    double sqr(double x) { return x*x; }
    double dist2(double x1, double y1, double x2, double y2) { return sqr(x1-x2)+sqr(y1-y2); }
    double dist(double x1, double y1, double x2, double y2) { return sqrt(dist2(x1,y1,x2,y2)); }
    double len2(double x, double y) { return sqr(x)+sqr(y); }
    double len(double x, double y) { return sqrt(len2(x,y)); }



    double SolveTomaz(int X1, int Y1, int X2, int Y2, int xc, int yc, int r) {

        auto ScoreTomaz = [=] (double x, double y) -> double {
            double vx=x-xc, vy=y-yc;
            double l=len(vx,vy);
            double xp=xc+r*vx/l, yp=yc+r*vy/l;
            return dist(X1,Y1,xp,yp)+dist(xp,yp,X2,Y2);
        };

        int d1=dist2(X1,Y1,xc,yc), d2=dist2(X2,Y2,xc,yc);
        int r2=r*r;
        int in1=d1<r2, in2=d2<r2;
        // touching the circle or inside/outside?
        if (d1==r2 || d2==r2 || in1!=in2) {
            return dist(X1,Y1,X2,Y2);
        }
        // straight line intersects the circle?
        double ux=xc-X1, uy=yc-Y1;
        double vx=X2-X1, vy=Y2-Y1;
        double proj=(ux*vx+uy*vy)/len2(vx,vy);
        double px=X1+proj*vx, py=Y1+proj*vy;
        if (!in1 && !in2 && 0<=proj && proj<=1 && dist2(xc,yc,px,py)<=r2) {
            return dist(X1,Y1,X2,Y2);
        }
        // ternary search along the line (X1,Y1) - (X2,Y2)
        auto Angle = [=] (double lambda) -> double { return atan2(Y1 + lambda * vy - yc, X1 + lambda * vx - xc); };
        auto AScore = [=] (double lambda) -> double { return ScoreTomaz(X1 + lambda * vx, Y1 + lambda * vy); };
        double a=0, b=1;
        for (int it=0;it<100;it++) {
            double k1=(2*a+b)/3, k2=(a+2*b)/3;
            if (false) printf("Tomaz it %d] %.9f %.9f %.9f %.9f, values %.9f %.9f %.9f %.9f\n", it, Angle(a), Angle(k1), Angle(k2), Angle(b), AScore(a), AScore(k1), AScore(k2), AScore(b));
            if (ScoreTomaz(X1+k1*vx, Y1+k1*vy) > ScoreTomaz(X1+k2*vx, Y1+k2*vy)) a=k1;
            else b=k2;
        }
        return ScoreTomaz(X1+a*vx, Y1+a*vy);
    }

}

struct Case
{
    int caseNo;
    int xA, xB, xC, yA, yB, yC, r;
    void Read()
    {
        int ok_ = scanf("%d %d %d %d %d %d %d", &xA, &yA, &xB, &yB, &xC, &yC, &r);
        Assert(ok_ == 7);
        Assert(-MaxCoord <= xA); Assert(xA <= MaxCoord);
        Assert(-MaxCoord <= yA); Assert(yA <= MaxCoord);
        Assert(-MaxCoord <= xB); Assert(xB <= MaxCoord);
        Assert(-MaxCoord <= yB); Assert(yB <= MaxCoord);
        Assert(-MaxCoord <= xC); Assert(xC <= MaxCoord);
        Assert(-MaxCoord <= yC); Assert(yC <= MaxCoord);
        Assert(0 <= r); Assert(r <= MaxR);
    }

    template<typename T>
    static inline long double SqrtDist(T dx, T dy) { return sqrt((long double) (dx * dx + dy * dy)); }
    template<typename T>
    static inline long double SqrtDist(T x1, T y1, T x2, T y2) { return SqrtDist(x1 - x2, y1 - y2); }

    double TernarySearch()
    {
        if (r == 0) return SqrtDist(xA, yA, xC, yC) + SqrtDist(xB, yB, xC, yC);
        if (xA == xB && yA == yB) return 2 * abs(SqrtDist(xA, yA, xC, yC) - r);
        if (xA == xC && yA == yC || xB == xC && yB == yC) {
            double d = SqrtDist(xA, yA, xB, yB); return d >= r ? d : 2 * r - d; }
        // Move the coordinate origin into the centre of the circle.
        long double xA_ = xA - xC, yA_ = yA - yC, xB_ = xB - xC, yB_ = yB - yC;
        long double angleA = atan2(yA_, xA_), angleB = atan2(yB_, xB_);
        if (angleB < angleA - pi) angleB += 2 * pi; 
        else if (angleB > angleA + pi) angleB -= 2 * pi;
        if (angleB < angleA) swap(angleA, angleB);
        if (angleA < 0) { angleB += 2 * pi; angleA += 2 * pi; }
        auto Eval = [this, xA_, yA_, xB_, yB_] (long double phi) { long double x = r * cos(phi), y = r * sin(phi); return SqrtDist(xA_, yA_, x, y) + SqrtDist(xB_, yB_, x, y); };
        long double valueA = Eval(angleA), valueB = Eval(angleB);
        for (int i = 0; i < 100; ++i)
        {
            long double angleM1 = (2 * angleA + angleB) / 3.0, angleM2 = (angleA + 2 * angleB) / 3.0;
            long double valueM1 = Eval(angleM1), valueM2 = Eval(angleM2);
            if (false) fprintf(stderr, "%d] Ternary %.12Lf %.12Lf %.12Lf %.12Lf [delta %Lg] values %.9Lf %.9Lf %.9Lf %.9Lf\n", i, angleA, angleM1, angleM2, angleB, angleB - angleA, valueA, valueM1, valueM2, valueB);
            long double dOld = angleB - angleA;
            if (valueM1 > valueM2) angleA = angleM1, valueA = valueM1;
            else angleB = angleM2, valueB = valueM2;
            long double dNew = angleB - angleA;
            if (dNew >= dOld) break;
        }
        long double angle = (angleA + angleB) / 2.0;
        long double xD = r * cos(angle), yD = r * sin(angle);
        long double result = SqrtDist(xA_, yA_, xD, yD) + SqrtDist(xB_, yB_, xD, yD);
        if (false) printf("direct: phi = %.12Lf, f = %.12Lf\n", angle, result);
        return result;
    }

    double VidBrutal()
    {
        long double xA_ = xA - xC, yA_ = yA - yC, xB_ = xB - xC, yB_ = yB - yC;
        auto Eval = [this, xA_, yA_, xB_, yB_] (long double phi) { long double x = r * cos(phi), y = r * sin(phi); return SqrtDist(xA_, yA_, x, y) + SqrtDist(xB_, yB_, x, y); };
        long double lb = 0, ub = 2 * pi;
        for (int i = 0; i < 50; ++i)
        {
            int n = (i == 0) ? 5000 : 10;
            long double margin = (ub - lb) / ((long double) n);
            long double fBest = 0, phiBest = 0;
            for (int j = 0; j < n; ++j)
            {
                long double phi = lb + margin * j;
                long double f = Eval(phi);
                if (j == 0 || f < fBest) phiBest = phi, fBest = f;
            }
            if (false) printf("i = %d: lb = %.12Lf, ub = %.12Lf [delta %Lg], best phi = %.9Lf, f = %.9Lf\n", i, lb, ub, ub - lb, phiBest, fBest);
            lb = phiBest - margin; ub = phiBest + margin;
        }
        return Eval(lb);
    }

    double GradientDescent(double &phi)
    {
        double coef = 0.1, f = -1, df, fPrev = -1;
        int nSteps = 0, nSinceReset = 0;
        if (false) printf("A(%d, %d), B(%d, %d), C(%d, %d), r = %d\n", xA, yA, xB, yB, xC, yC, r);
        for (nSteps = 0; nSteps <= 10000; ++nSteps, ++nSinceReset)
        {
            // f(phi) = sum_i sqrt( (xi - xC - r cos phi)^2 + (yi - yC - r sin phi)^2)
            // f'(phi) = sum_i (1/(2 sqrt(...))) [2 (xi - xC - r cos phi) r sin phi + 2 (yi - yC - r sin phi) (-r cos phi)]
            //         = sum_i (r/sqrt(...)) [(xi - xC - r cos phi) sin phi - (yi - yC - r sin phi) cos phi]
            f = 0; df = 0;
            double f1_ = 0, f2_ = 0;
            for (int i = 0; i < 2; ++i)
            {
                double xi = (i == 0) ? xA : xB, yi = (i == 0) ? yA : yB;
                double Cos = cos(phi), Sin = sin(phi);
                double xPart = xi - xC - r * Cos, yPart = yi - yC - r * Sin;
                double Sqrt = sqrt(xPart * xPart + yPart * yPart);
                if (i == 0) f1_ = Sqrt; else f2_ = Sqrt;
                f += Sqrt;
                df += (r / max(Sqrt, 1e-6)) * (xPart * Sin - yPart * Cos);
            }
            /*if (coef > 1e-6) coef *= 0.99; */
            //if (nSteps % 30 == 0) coef *= 0.2;
            if (nSteps > 0 && f > fPrev && nSinceReset > 10) coef *= 0.2;
            //double Coef = coef; if (abs(df) > 1e-8) Coef = min(Coef, abs(f) / abs(df) / 2.0);
            //phi += (df * coef;
            if (verbose) if (nSteps % 1 == 0) printf("Step %d: phi = %.6f, step = %.6g, f = %.6f + %.6f = %.6f [%g], df = %.6g\n", nSteps, phi, coef, f1_, f2_, f, abs(f - fPrev), df);
            phi += (df > 0 ? -1 : 1) * coef;
            while (phi < -pi) phi += 2 * pi;
            while (phi > pi) phi -= 2 * pi;
            if (abs(df) < 1e-10) break;
            if (nSteps > 0 && abs(f - fPrev) < 1e-10 * f) break;
            if (nSteps > 0 && abs(f - fPrev) < 1e-10) break;
            fPrev = f;
        }
        if (verbose) printf("%d steps: phi = %.8f, step = %.8g, f = %.8f, df = %.8f\n", nSteps, phi, coef, f, df);
        return f;
    }

    double GradientDescent()
    {
        constexpr int nTries = 3;
        double fMin = -1, fMax = -1;
        for (int tryNo = 0; tryNo < nTries; ++tryNo)
        {
            double phi0 = tryNo / double(nTries) * (2 * pi) + sqrt(2);
            double phi = phi0; double f = GradientDescent(phi);
            if (fMin < 0 || f < fMin) fMin = f;
            if (fMax < 0 || f > fMax) fMax = f;
            if (verbose) printf("Try %d: phi %.6f -> %.6f, f = %.6f\n", tryNo, phi0, phi, f);
        }
        if (false) if (fMax - fMin > 1e-6) fprintf(stderr, "case %d: fMax - fMin = %g\n", caseNo, fMax - fMin);
        double f2 = TernarySearch();
        double f3 = VidBrutal();
        double f4 = Tomaz::SolveTomaz(xA, yA, xB, yB, xC, yC, r);
        if (abs(f2 - fMin) > 1e-6 || abs(f2 - f3) > 1e-6 || abs(f3 - f4) > 1e-6) fprintf(stderr, "case %d: grad desc %.9g, ternary %.9g, Vid %.9g [delta = %.9g], Tomaz %.9f [delta %.9g]\n", caseNo, fMin, f2, f3, f3 - f2, f4, f3 - f4);
        return fMin;
    }

    template<typename T>
    static inline bool InRange(T t1, T t, T t2) 
    {
        while (t1 < 0) t1 += 2 * pi; while (t1 >= 2 * pi) t1 -= 2 * pi;
        while (t2 < 0) t2 += 2 * pi; while (t2 >= 2 * pi) t2 -= 2 * pi;
        while (t < 0) t += 2 * pi; while (t >= 2 * pi) t -= 2 * pi;
        if (t1 <= t2) return t1 <= t && t <= t2; 
        else return (t1 <= t) || (t <= t2);
    }

    template<typename T>
    static inline void DistToBbox(T x1, T y1, T x2, T y2, T x, T y, T& minDist, T& maxDist)
    {
        if (x1 > x2) swap(x1, x2);
        if (y1 > y2) swap(y1, y2);
        T dx = (x < x1) ? (x1 - x) : (x > x2) ? (x - x2) : 0;
        T dy = (y < y1) ? (y1 - y) : (y > y2) ? (y - y2) : 0; 
        minDist = sqrt(dx * dx + dy * dy);
        dx = max(abs(x - x1), abs(x - x2)); dy = max(abs(y - y1), abs(y - y2));
        maxDist = sqrt(dx * dx + dy * dy);
        if (minDist > maxDist) swap(minDist, maxDist); // this can happen due to numerical inaccuracies when (x1, y1) and (x2, y2) are extremely close together
    }

    // This assumes that the bounds on the denominator are positive.
    template<typename T>
    static inline T lbRatio(T num, T lbDen, T ubDen) { Assert(0 < lbDen); Assert(lbDen <= ubDen); return num < 0 ? num / lbDen : num / ubDen; }
    template<typename T>
    static inline T ubRatio(T num, T lbDen, T ubDen) { Assert(0 < lbDen); Assert(lbDen <= ubDen); return num < 0 ? num / ubDen : num / lbDen; }

    double VidVariant() const
    {
        //fprintf(stderr, "%d %d %d %d %d %d %d\n", xA, yA, xB, yB, xC, yC, r); fflush(stderr);
        long double xA_ = xA - xC, yA_ = yA - yC, xB_ = xB - xC, yB_ = yB - yC;
        long double dA = SqrtDist(xA_, yA_), dB = SqrtDist(xB_, yB_);
        long double uA = atan2(yA_, xA_), uB = atan2(yB_, xB_);
        if (r == 0) return dA + dB;
        if (dA <= r && dB >= r || dA >= r && dB <= r) return SqrtDist(xA_, yA_, xB_, yB_);
        xA_ /= r; xB_ /= r; yA_ /= r; yB_ /= r; dA /= r; dB /= r;
        //
        struct Cand { long double t1, t2, f, df, lbF; };
        vector<Cand> cands; long double fMin = max((long double) 2, dA + dB); // 2 is useful if they are both inside the circle; dA + dB if the are both outside
        int totalCands = 0, maxCandsAtATime = 0;
        enum { verbose = false };
        auto AddCand = [=, &cands, &fMin, &totalCands] (long double t1, long double t2) 
        {
            long double t = (t1 + t2) / 2;
            long double xT = cos(t), yT = sin(t);
            long double dTA = SqrtDist(xT, yT, xA_, yA_), dTB = SqrtDist(xT, yT, xB_, yB_);
            long double f = dTA + dTB; ++totalCands;
            // f(t) = sum_i sqrt[(cos t - x_i)^2 + (sin t - y_i)^2]
            // f'(t) = sum_i 1/sqrt[(cos t - x_i)^2 + (sin t - y_i)^2] * ((cos t - x_i) (-sin t) + (sin t - y_i) cos t)
            //       = sum_i [x_i sin t - y_i cos t] / sqrt[(cos t - x_i)^2 + (sin t - y_i)^2]
            // Let us represent (x_i, y_i) in polar coordinates as x_i = z_i cos u_i and y_i = z_i sin u_i.
            // f'(t) = sum_i z_i [cos u_i sin t - sin u_i cos t] / sqrt[(cos t - x_i)^2 + (sin t - y_i)^2]
            //       = sum_i z_i sin(t - u_i) / sqrt[(cos t - x_i)^2 + (sin t - y_i)^2]
            // For bounds on f'(t): 
            // - The term sin(t - u_i) achieves its minimum of -1 at t = u_i - pi/2 and its maximum of 1 at u_i + pi/2.
            //   So the minimum and maximum value of this term for t in [t1, t2] will be either at t1 or at t2 
            //   or at u_i +/- pi/2, if this value lies in the range [t1, t2].
            // - For the term sqrt[...], we can get a lower and an upper bound using the bounding box
            //   of the points (cos t1, sin t1) and (cos t2, sin t2); and we can also get a lower and an upper
            //   bound by using the distance between A (or B) and the closest/farthest point on the unit circle.
            long double lbSinA = min(sin(t1 - uA), sin(t2 - uA)); if (InRange(t1, uA - pi * 0.5, t2)) lbSinA = min(lbSinA, (long double) -1);
            long double lbSinB = min(sin(t1 - uB), sin(t2 - uB)); if (InRange(t1, uB - pi * 0.5, t2)) lbSinB = min(lbSinB, (long double) -1);
            long double ubSinA = max(sin(t1 - uA), sin(t2 - uA)); if (InRange(t1, uA + pi * 0.5, t2)) ubSinA = max(lbSinA, (long double) 1);
            long double ubSinB = max(sin(t1 - uB), sin(t2 - uB)); if (InRange(t1, uB + pi * 0.5, t2)) ubSinB = max(lbSinB, (long double) 1);
            long double x1 = cos(t1), y1 = sin(t1), x2 = cos(t2), y2 = sin(t2);
            long double f1 = SqrtDist(x1, y1, xA_, yA_) + SqrtDist(x1, y1, xB_, yB_);
            long double f2 = SqrtDist(x2, y2, xA_, yA_) + SqrtDist(x2, y2, xB_, yB_);
            long double lbSqrtA, ubSqrtA; DistToBbox(x1, y1, x2, y2, xA_, yA_, lbSqrtA, ubSqrtA);
            long double lbSqrtB, ubSqrtB; DistToBbox(x1, y1, x2, y2, xB_, yB_, lbSqrtB, ubSqrtB);
            lbSqrtA = max(lbSqrtA, min(abs(dA - 1), abs(dA + 1))); ubSqrtA = min(ubSqrtA, max(abs(dA - 1), abs(dA + 1)));
            lbSqrtB = max(lbSqrtB, min(abs(dB - 1), abs(dB + 1))); ubSqrtB = min(ubSqrtB, max(abs(dB - 1), abs(dB + 1)));
            if (lbSqrtA > ubSqrtA) swap(lbSqrtA, ubSqrtA); // this can happen due to numerical inaccuracies
            if (lbSqrtB > ubSqrtB) swap(lbSqrtB, ubSqrtB); // when the range [t1, t2] is very narrow
            long double lbDF = dA * lbRatio(lbSinA, lbSqrtA, ubSqrtA) + dB * lbRatio(lbSinB, lbSqrtB, ubSqrtB);
            long double ubDF = dA * ubRatio(ubSinA, lbSqrtA, ubSqrtA) + dB * ubRatio(ubSinB, lbSqrtB, ubSqrtB);
            Assert(lbDF <= ubDF);
            // We use the lower bound on f' for the range [t, t2] and the upper bound for the range [t1, t].
            // Thus we could have made the bounds tighter by using only these ranges instead of 
            // calculating both bounds for the wider range [t1, t2].
            long double dt = (t2 - t1) / 0.5;
            long double lbF = f + min((long double) 0, min(dt * lbDF, - dt * ubDF));
            lbF = max(lbF, f1 + 2 * dt * lbDF); lbF = max(lbF, f2 - 2 * dt * ubDF);
            if (verbose) fprintf(stderr, "Cand: t = [%.9Lf %.9Lf] %.9Lf -> f = %.9Lf, bounds on df %.9Lg %.9Lg, lbound on f = %.9Lf [delta %.9Lg]\n", t1, t2, t, f, lbDF, ubDF, lbF, lbF - f);
            if (lbF > fMin - 1e-8 / r) return; // the interval [t1, t2] is hopeless
            if (f < fMin) fMin = f;
            long double df = (xA * yT - yA * xT) / dTA + (xB * yT - yB * xT) / dTB;
            cands.push_back({t1, t2, f, df, lbF}); 
        };
        int nSteps = 10, fanout = 5;
        for (int i = 0; i < nSteps; ++i) AddCand((i * 2 * pi) / ((long double) nSteps), ((i + 1) * 2 * pi) / ((long double) nSteps));
        vector<Cand> oldCands; long double stepSize = 2 * pi / nSteps;
        int nIter = 0;
        for ( ; stepSize > 1e-16 && ! cands.empty(); ++nIter)
        {
            stepSize /= fanout;
            maxCandsAtATime = max(maxCandsAtATime, (int) cands.size());
            swap(oldCands, cands); cands.clear();
            long double lbMin = fMin; for (auto & cand : oldCands) lbMin = min(lbMin, cand.lbF);
            if (verbose) fprintf(stderr, "Iter %d: %d cands, fMin = %.9Lf, lbF = %.9Lf (delta %.9Lg)\n", nIter, int(oldCands.size()), fMin, lbMin, lbMin - fMin);
            for (auto &cand : oldCands)
            {
                if (cand.lbF > fMin) continue;
                /*
                long double tMid = (cand.t1 + cand.t2) * 0.5;
                AddCand(cand.t1, tMid); AddCand(tMid, cand.t2);
                */
               for (int i = 0; i < fanout; ++i) AddCand((cand.t1 * (fanout - i) + cand.t2 * i) / ((long double) fanout), (cand.t1 * (fanout - (i + 1)) + cand.t2 * (i + 1)) / ((long double) fanout));
            }
            lbMin = fMin; for (auto & cand : cands) lbMin = min(lbMin, cand.lbF);
            if (verbose) fprintf(stderr, "Iter %d: %d -> %d cands, stepSize = %.9Lg, fMin = %.9Lf, lbF = %.9Lf (delta %.9Lg)\n", nIter, int(oldCands.size()), int(cands.size()), stepSize, fMin, lbMin, lbMin - fMin);
            if (lbMin >= fMin - 1e-8 / r) break; // further improvements are to small to be worth bothering with
        }
        fprintf(stderr, "%d iterations, %d cands (max %d at a time)\n", nIter, totalCands, maxCandsAtATime);
        return fMin * r;
    }

    template<typename T>
    static inline T LinInterp(T x1, T y1, T x2, T y2, T x) { return (x1 == x2) ? (y1 + y2) * 0.5 : ((x2 - x) * y1 + (x - x1) * y2) / (x2 - x1); }

    double NewtonStyle() const
    {
        long double xA_ = xA - xC, yA_ = yA - yC, xB_ = xB - xC, yB_ = yB - yC;
        long double dA = SqrtDist(xA_, yA_), dB = SqrtDist(xB_, yB_);
        long double uA = atan2(yA_, xA_), uB = atan2(yB_, xB_);
        if (r == 0) return dA + dB;
        if (dA <= r && dB >= r || dA >= r && dB <= r) return SqrtDist(xA_, yA_, xB_, yB_);
        if (! MUST_TOUCH_EDGE) {
            if (dA <= r && dB <= r) return SqrtDist(xA_, yA_, xB_, yB_);
            long double xAB = xB_ - xA_, yAB = yB_ - yA_;
            long double lambda = (xA_ * xAB + yA_ * yAB) / (xAB * xAB + yAB * yAB);  if (lambda < 0) lambda = 0; else if (lambda > 1) lambda = 1;
            long double xProj = xA_ * (1 - lambda) + xB_ * lambda, yProj = yA_ * (1 - lambda) + yB_ * lambda;
            if (xProj * xProj + yProj * yProj <= r * r) return SqrtDist(xA_, yA_, xB_, yB_);
        }
        xA_ /= r; xB_ /= r; yA_ /= r; yB_ /= r; dA /= r; dB /= r;
        auto Calc = [=] (long double t, long double &f, long double &df) {
            long double x = cos(t), y = sin(t);
            long double dTA = SqrtDist(x, y, xA_, yA_), dTB = SqrtDist(x, y, xB_, yB_);
            f = dTA + dTB; df = (xA_ * y - yA_ * x) / dTA + (xB_ * y - yB_ * x) / dTB; };
        long double fMin = SqrtDist(xA_, yA_, (long double) 1, (long double) 0) + SqrtDist(xB_, yB_, (long double) 1, (long double) 0); // any point on the circle will do
        long double tPrev = -1, fPrev = -1, dfPrev = -1;
        constexpr int n = 100;
        for (int i = 0; i <= n; ++i)
        {
            long double t = 2 * pi * i / ((long double) n);
            long double f, df; Calc(t, f, df);
            if (i > 0 && (dfPrev <= 0 && 0 <= df || dfPrev >= 0 && 0 >= df))
            {
                long double t1 = tPrev, f1 = fPrev, df1 = dfPrev, t2 = t, f2 = f, df2 = df;
                //for (int steps = 0; steps < 100; ++steps)
                while (true)
                {
                    //long double te = LinInterp(df1, t1, df2, t2, (long double) 0);  // this one has problems due to inaccuracy on '1 1000 2 1000 0 0 1000'
                    long double te = (t1 + t2) / 2;
                    long double fe, dfe; Calc(te, fe, dfe);
                    if (false) fprintf(stderr, "%.9Lf %.9Lf %+.2Le | %.9Lf %.9Lf %+.2Le | %.9Lf %.9Lf %+.2Le\n", t1, f1, df1, te, fe, dfe, t2, f2, df2);
                    if (fe < fMin) fMin = fe;
                    if (abs(dfe) < 1e-15 || abs(t2 - t1) < 1e-15) break;
                    if (df1 <= 0 && 0 <= dfe || df1 >= 0 && 0 >= df1) t2 = te, f2 = fe, df2 = dfe;
                    else { Assert(dfe <= 0 && 0 <= df2 || dfe >= 0 && 0 >= df2); t1 = te, f1 = fe, df1 = dfe; }
                }
            }
            tPrev = t; fPrev = f; dfPrev = df;
        }
        return fMin * r;
    }

};

void Test(int seed = 0)
{
    double maxDelta = 0, maxDelta2 = 0;
    for (long long int nIter = 0; ; ++nIter)
    {
        if (nIter % 1000 == 0) { fprintf(stderr, "%lld  %.9g %.g      \r", nIter, maxDelta, maxDelta2); fflush(stderr); }
        mt19937_64 r(123 + nIter + 1'000'000'000LL * seed);
        Case c; c.caseNo = (int) nIter;
        c.xC = 0; c.yC = 0; c.r = 30;
        c.xA = uniform_int_distribution(-200, 200)(r);
        c.xB = uniform_int_distribution(-200, 200)(r);
        c.yA = uniform_int_distribution(-200, 200)(r);
        c.yB = uniform_int_distribution(-200, 200)(r);
        double f2 = c.TernarySearch(), f3 = c.VidBrutal(), f4 = Tomaz::SolveTomaz(c.xA, c.yA, c.xB, c.yB, c.xC, c.yC, c.r);
        double fTested = f2;
        double d = abs(fTested - f3);
        if ((fTested > f3 && fTested - f3 > maxDelta) || (f3 > fTested && f3 - fTested > maxDelta2))
        {
            printf("case %lld: A(%d, %d), B(%d, %d), C(%d, %d), r = %d:  ternary %.12g, Tomaz %.12g, Vid %.12g, abs = %g, rel = %g\n", nIter,
                c.xA, c.yA, c.xB, c.yB, c.xC, c.yC, c.r, f2, f4, f3, fTested - f3, (fTested > f3 ? fTested / f3 - 1 : f3 / fTested - 1));
            fflush(stdout);
            if (fTested > f3 && abs(fTested / f3 - 1) < 1e-3) maxDelta = d;
            if (f3 > fTested) maxDelta2 = d;
        }


    }
}

void Test2()
{
    long long int caseNo = 0; double maxDelta = 0;
    for (int u = 1; u <= 1000; ++u)
        for (int xA = -u; xA <= u; ++xA) for (int yA = -u; yA <= u; ++yA)
        {
            if (xA != u && xA != -u && yA != u && yA != -u) continue;
            for (int xB = -u; xB <= u; ++xB) for (int yB = -u; yB <= u; ++yB)
            {
                Case c; c.caseNo = (int) caseNo; 
                if (++caseNo % 1000 == 0) { fprintf(stderr, "%lld: u = %d, A(%d, %d), B(%d, %d)            \r", caseNo, u, xA, yA, xB, yB); fflush(stderr); }
                c.xC = 0; c.yC = 0; c.r = 30;
                c.xA = xA; c.xB = xB; c.yA = yA; c.yB = yB;
                double f2 = c.TernarySearch(), f3 = c.VidBrutal();
                double d = abs(f2 - f3);
                if (d > maxDelta)
                {
                    printf("case %lld: A(%d, %d), B(%d, %d), C(%d, %d), r = %d:  ternary %.12g, Vid %.12g, abs = %g, rel = %g\n", caseNo,
                        c.xA, c.yA, c.xB, c.yB, c.xC, c.yC, c.r, f2, f3, f2 - f3, f2/f3 - 1);
                    fflush(stdout);
                    maxDelta = d;
                }
            }

        }
}

int main(int argc, char** argv)
{
    //Test(argc > 1 ? atoi(argv[1]) : 0);
    //Test2();
    //freopen("moon97.in", "rt", stdin);
    int T; 
    int ok_ = scanf("%d", &T); Assert(ok_ == 1); Assert(1 <= T); Assert(T <= MaxT);
    for (int caseNo = 0; caseNo < T; ++caseNo)
    {
        Case c; c.caseNo = caseNo; c.Read();
        //printf("%.15g\n", c.GradientDescent());
        //printf("%.15g\n", c.VidVariant());
        printf("%.15g\n", c.NewtonStyle());
        //
    }
    return 0;
}