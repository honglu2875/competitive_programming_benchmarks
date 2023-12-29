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
#define DEBUG_ITER(x) cerr << #x << " = ["; for (auto _ : x) cerr << _ << "| "; cerr << "]" << endl;
template<class A, class B> ostream& operator<<(ostream& os, pair<A,B> &p) {
	os << "(" << p.first << ", " << p.second << ")";
	return os;
}

typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int,int> PII;
typedef vector<int> VI;
typedef array<int,3> III;
typedef vector<string> VS;

#define N 600
#define D 1000000

int n,m,d,start;
string s;

int64 h[D+1], pf[D+1], pfi[D+1];
int64 hf=31, hfi=129032259, MOD=1000000007;

typedef PII Sub;
typedef vector<PII> Chain;

int sign(int x) {
	if (x<0) return -1;
	else if (x>0) return 1;
	else return 0;
}

inline uint64 subhash(int i, int n) {
	return ((h[i+n-1]-h[i-1]+MOD)*pfi[i])%MOD;
}

int cmp(int ai, int an, int bi, int bn) {
	int n=min(an,bn);
	int lo=0, hi=n+1;
	if (subhash(ai,n)==subhash(bi,n)) lo=n;
	while (lo+1<hi) {
		int mid=(lo+hi)/2;
		if (subhash(ai,mid)==subhash(bi,mid)) lo=mid;
		else hi=mid;
	}
	if (lo<n) return sign(s[ai+lo]-s[bi+lo]);
	else {
		if (an==bn) return 0;
		else if (an<bn) return -1;
		else return 1;
	}
}

int cmp(const Chain &a, const Chain &b) {
	int i=(int)a.size()-1,j=(int)b.size()-1,ai,bi,an=0,bn=0;
	while (1) {
		while (i>=0 && an==0) { ai=a[i].first; an=a[i].second; i--; }
		while (j>=0 && bn==0) { bi=b[j].first; bn=b[j].second; j--; }
		if (an==0 && bn==0) return 0;
		else if (an==0) return -1;
		else if (bn==0) return 1;
		int m=min(an,bn);
		int c=cmp(ai,m,bi,m);
		if (c!=0) return c;
		ai+=m; an-=m; bi+=m; bn-=m;
	}
}

bool operator<(const Chain &a, const Chain &b) { return cmp(a,b)<0; }
bool operator!=(const Chain &a, const Chain &b) { return cmp(a,b)!=0; }

int outdeg[N+1];
vector<int> prv[N+1];
vector<pair<int,Sub>> adj[N+1];

int done[N+1];
Chain shortest[N+1];
int nxt[N+1];
vector<int> paths[N+1];

int main() {
	//freopen("test.in","r",stdin);
	cin >> n >> m >> d >> start;
	cin >> s;
	s="#"+s;

	h[0]=s[0]-'a';
	pf[0]=1; pfi[0]=1;
	FOR (i,1,d) {
		pf[i]=(pf[i-1]*hf)%MOD; pfi[i]=(pfi[i-1]*hfi)%MOD;
		h[i]=(h[i-1]+pf[i]*(s[i]-'a'))%MOD;
	}

	REP (i,m) {
		int u,v,p,l;
		cin >> u >> v >> p >> l;
		Sub e = {p,l};
		adj[u].push_back({v,e});
		prv[v].push_back(u);
		outdeg[u]++;
	}

	vector<int> topo;
	queue<int> term;
	FOR (i,1,n) if (outdeg[i]==0) term.push(i);
	while (!term.empty()) {
		int x=term.front(); term.pop();
		topo.push_back(x);
		for (int y : prv[x]) {
			outdeg[y]--;
			if (outdeg[y]==0) term.push(y);
		}
	}

	FOR (fin,1,n) {
		// check reachability
		INIT(done,0);
		done[fin]=1;
		for (int x : topo) {
			for (auto& [y,e] : adj[x]) if (done[y]) done[x]=1;
		}
		if (!done[start]) continue;
		// compute shortest paths
		INIT(done,0); INIT(nxt,0);
		done[fin]=1; shortest[fin]={};
		for (int x : topo) {
			for (auto [y,e] : adj[x]) if (done[y]) {
				Chain &c=shortest[y];
				c.push_back(e);
				if (!done[x] || c<shortest[x]) {
					done[x]=1;
					shortest[x]=c;
					nxt[x]=y;
				}
				c.pop_back();
			}
		}
		// reconstruct
		vector<int> &path = paths[fin];
		int x=start;
		while (x!=0) {
			path.push_back(x);
			x=nxt[x];
		}
	}

	FOR (i,1,n) {
		vector<int> &path = paths[i];
		printf("%d",path.size());
		for (int x : path) printf(" %d",x);
		printf("\n");
	}
	return 0;
}
