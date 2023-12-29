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

typedef pair<int,int> PII;
typedef long long int64;
typedef vector<int> VI;
typedef vector<int64> VL;
typedef vector<PII> VII;
typedef vector<double> VD;

#define N 500

int n;
int64 cap[N][N];
int source, sink;
int64 inf=1LL<<60;

int vis[N], lim, mark;
int64 flow(int x, int64 f) {
    if (x==sink) return f;
    vis[x]=mark;
    REP (y,n) {
        if (cap[x][y]>=lim && vis[y]!=mark) {
            int64 ff=flow(y,min(f,cap[x][y]));
            if (ff!=0) {
                cap[x][y]-=ff;
                cap[y][x]+=ff;
                return ff;
            }
        }
    }
    return 0;
}

int main() {
	int I,C,T;
	cin >> I >> C >> T;
	VI c(C), x(I), t(T);
	REP (i,C) cin >> c[i];
	REP (i,I) cin >> x[i];
	REP (i,T) cin >> t[i];
	VL ing(C);
	REP (i,C) {
		REP (j,I) {
			int a;
			cin >> a;
			ing[i]+=(int64)a*x[j];
		}
	}
	vector<VI> req(C);
	REP (i,C) {
		int k;
		cin >> k;
		REP (j,k) {
			int b;
			cin >> b;
			req[i].PB(b-1);
		}
	}
	source=C+T; sink=source+1; n=C+T+2;
	int64 total=0;
	REP (i,C) {
		int64 profit=max(c[i]-ing[i], 0LL);
		total+=profit;
		cap[source][i]=profit;
	}
	REP (i,T) {
		cap[C+i][sink]=t[i];
	}
	REP (i,C) {
		for (int j : req[i]) {
			cap[i][C+j]=inf;
		}
	}
	int64 f=0;
	for (lim=1<<30;lim>=1;lim/=2) {
		while (1) {
			mark++;
			int64 df=flow(source, inf);
			if (df==0) break;
			f+=df;
		}
	}
	cout << total - f << endl;
	return 0;
}
