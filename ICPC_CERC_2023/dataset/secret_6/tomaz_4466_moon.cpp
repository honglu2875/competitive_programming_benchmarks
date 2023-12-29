#include <bits/stdc++.h>
using namespace std;

#define ALL(c) (c).begin(),(c).end()
#define PB push_back
#define IN(x,c) (find(c.begin(),c.end(),x) != (c).end())
#define REP(i,n) for (int i=0;i<(int)(n);i++)
#define FOR(i,a,b) for (int i=(a);i<=(b);i++)
#define INIT(a,v) memset(a,v,sizeof(a))
template<class A, class B> A cvt(B x) { stringstream ss; ss<<x; A y; ss>>y; return y; }

#define SPC << " " <<
#define DEBUG(x) cerr << #x << " = "; cerr << x << endl;
#define DEBUG_ITER(x) cerr << #x << " = "; for (auto _ : x) cerr << _ << ' '; cerr << endl;

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

int tests;
int X1,Y1,X2,Y2,xc,yc,r;

double score(double x, double y) {
	double vx=x-xc, vy=y-yc;
	double l=len(vx,vy);
	double xp=xc+r*vx/l, yp=yc+r*vy/l;
	return dist(X1,Y1,xp,yp)+dist(xp,yp,X2,Y2);
}

int main() {
	//freopen("circle.in","r",stdin);
	//freopen("circle.out","w",stdout);
	cin >> tests;
	REP (test,tests) {
		cin >> X1 >> Y1 >> X2 >> Y2 >> xc >> yc >> r;
		int d1=dist2(X1,Y1,xc,yc), d2=dist2(X2,Y2,xc,yc);
		int r2=r*r;
		int in1=d1<=r2, in2=d2<=r2;
		// already inside the circle
		if (in1 || in2) {
			printf("%.9lf\n",dist(X1,Y1,X2,Y2));
			continue;
		}
		// straight line intersects the circle?
		double ux=xc-X1, uy=yc-Y1;
		double vx=X2-X1, vy=Y2-Y1;
		double proj=(ux*vx+uy*vy)/len2(vx,vy);
		double px=X1+proj*vx, py=Y1+proj*vy;
		if (0<=proj && proj<=1 && dist2(xc,yc,px,py)<=r2) {
			printf("%.9lf\n",dist(X1,Y1,X2,Y2));
			continue;
		}
		// ternary search along the line (X1,Y1) - (X2,Y2)
		double a=0, b=1;
		for (int it=0;it<100;it++) {
			double k1=(2*a+b)/3, k2=(a+2*b)/3;
			if (score(X1+k1*vx, Y1+k1*vy) > score(X1+k2*vx, Y1+k2*vy)) a=k1;
			else b=k2;
		}
		printf("%.9lf\n",score(X1+a*vx, Y1+a*vy));
	}
	return 0;
}
