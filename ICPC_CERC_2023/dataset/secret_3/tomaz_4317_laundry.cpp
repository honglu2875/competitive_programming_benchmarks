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

#define N 30000
#define L 300000
#define Q 300000

int n,q;
int d[N], fast[N], slow[N];
bitset<L+1> f[2];
int ans[Q];

int main() {
	scanf("%d %d",&n,&q);
	vector<array<int,3>> sheets;
	multiset<int> fast_times;
	fast_times.insert(0);
	int64 two=0, one=0;
	REP (i,n) {
		scanf("%d %d %d",&d[i],&fast[i],&slow[i]);
		sheets.push_back({slow[i],fast[i],d[i]});
		fast_times.insert(fast[i]);
		two+=d[i];
	}
	vector<PII> queries;
	REP (i,q) {
		int l;
		scanf("%d",&l);
		queries.push_back({l,i});
	}
	sort(ALL(queries)); reverse(ALL(queries));
	sort(ALL(sheets));
	INIT(ans,-1);
	int qi=0, t=*fast_times.rbegin();
	while (qi<(int)queries.size() && queries[qi].first>=two) {
		ans[queries[qi].second] = t; qi++;
	}
	f[0][L]=1;
	REP (i,n) {
		auto [slow_time, fast_time, x]=sheets[i];
		fast_times.erase(fast_times.find(fast_time));
		auto &row=f[i%2], &nxt=f[(i+1)%2];
		nxt = (row>>x) | row;
		two-=x; one+=x;
		int64 tail=one/2;
		if (tail>L) continue;
		auto cand = nxt>>(L-tail);
		int bi=cand._Find_first();
		if (bi>L) continue;
		int first=tail-bi, second=one-first;
		int t=max(slow_time, *fast_times.rbegin());
		while (qi<(int)queries.size() && queries[qi].first>=two+second) {
			ans[queries[qi].second] = t; qi++;
		}
	}
	REP (i,q) printf("%d\n",ans[i]);
	return 0;
}
