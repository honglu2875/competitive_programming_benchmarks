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

#define MAXN 600
#define MAXD 1000000

int n,m,d,start;
string s;

typedef PII Sub;
typedef vector<PII> Chain;

int lg2[MAXD+1];

namespace SA {
	typedef struct {
		int a; int b; int c;
	} triple;

	bool operator<(const triple &t1, const triple &t2) {
		if (t1.a==t2.a) {
			if (t1.b==t2.b) return t1.c<t2.c;
			else return t1.b<t2.b;
		} else return t1.a<t2.a;
	}

	int P[21][MAXD];  // 2^20 >= MAXD
	int N,k,sz;
	triple L[MAXD];
	int id2pos[MAXD], pos2id[MAXD];

	// stable count sort of L by first item for O(n*log n) construction
	int cnt[MAXD];
	triple L2[MAXD];
	void countSort(triple L[MAXD]) {
		for (int i=0;i<N;i++) cnt[i]=0;
		for (int i=0;i<N;i++) cnt[L[i].a]++;
		for (int i=1;i<N;i++) cnt[i]+=cnt[i-1];
		for (int i=N-1;i>=0;i--) L2[--cnt[L[i].a]] = L[i];
		for (int i=0;i<N;i++) L[i] = L2[i];
	}

	void suffixArray(string s) {
		N=s.size();
		for (int i=0;i<N;i++) P[0][i]=(int)s[i];
		for (k=1,sz=1; sz<N; k++,sz*=2) {
			for (int i=0;i<N;i++) { L[i].a=P[k-1][i]; L[i].b=(i+sz<N)?P[k-1][i+sz]:-1; L[i].c=i; }
			if (k==1) sort(L,L+N);  // not sorted by anything in first round
			else countSort(L);  // stable count sort by first item (already sorted by second item)
			for (int i=0;i<N;i++) {
				if (i!=0 && L[i].a==L[i-1].a && L[i].b==L[i-1].b) P[k][L[i].c]=P[k][L[i-1].c];
				else P[k][L[i].c]=i;
			}
		}
		for (int i=0;i<N;i++) {
			pos2id[i]=L[i].c;
			id2pos[L[i].c]=i;
		}
	}

	int lcp(int x, int y) {
		if (x==y) return N-x;
		int r=0;
		for (int len=k-1,size=1<<len; len>=0 && x<N && y<N; len--,size/=2) {
			if (P[len][x]==P[len][y]) { r+=size; x+=size; y+=size; }
		}
		return r;
	}


}

int sign(int x) {
	if (x<0) return -1;
	else if (x>0) return 1;
	else return 0;
}

int cmp(int ai, int an, int bi, int bn) {
	int n=min(an,bn);
	int l2=lg2[n], l=0;
	if (SA::P[l2][ai]==SA::P[l2][bi] && SA::P[l2][ai+n-(1<<l2)]==SA::P[l2][bi+n-(1<<l2)]) l=n;
	else l=SA::lcp(ai,bi);
	if (l<n) return sign(s[ai+l]-s[bi+l]);
	else return sign(an-bn);
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

int outdeg[MAXN+1];
vector<int> prv[MAXN+1];
vector<pair<int,Sub>> adj[MAXN+1];

int done[MAXN+1];
Chain shortest[MAXN+1];
int nxt[MAXN+1];
vector<int> paths[MAXN+1];

int main() {
	//freopen("test.in","r",stdin);
	//freopen("test.out","w",stdout);
	cin >> n >> m >> d >> start;
	cin >> s;
	SA::suffixArray(s);
	lg2[1]=0;
	FOR (i,2,d) lg2[i]=1+lg2[i/2];

	REP (i,m) {
		int u,v,p,l;
		cin >> u >> v >> p >> l;
		Sub e = {p-1,l};
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
		printf("%d",(int)path.size());
		for (int x : path) printf(" %d",x);
		printf("\n");
	}
	return 0;
}
