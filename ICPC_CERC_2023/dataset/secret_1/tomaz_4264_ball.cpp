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

double solve(vector<PII> p) {
	double d=0;
	assert(p.size()%2==0);
	REP (i,p.size()/2) {
		int j=(i+p.size()/2);
		d+=sqrt(pow(p[i].first-p[j].first,2)+pow(p[i].second-p[j].second,2));
	}
	return d;
}

int main() {
	int n;
	string s;
	cin >> n >> s;
	int nb=count(ALL(s),'B'), ng=count(ALL(s),'G');
	vector<PII> p(n);
	REP (i,n) cin >> p[i].first >> p[i].second;
	vector<PII> b, g;
	REP (i,n) {
		if (s[i]=='B') b.push_back(p[i]);
		if (s[i]=='G') g.push_back(p[i]);
	}
	printf("%.9lf\n",solve(b)+solve(g));
	return 0;
}
