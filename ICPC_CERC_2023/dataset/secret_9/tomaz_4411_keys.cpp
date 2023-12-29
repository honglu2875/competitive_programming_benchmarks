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

#define N 100000

int n,m;
vector<int> adj[N];
map<PII,int> edge;

int z=-1, az[N];
int vis[N],prv[N],depth[N];
void dfs(int x, int p=-1) {
	if (vis[x]) return;
	vis[x]=1; prv[x]=p;
	if (p!=-1) depth[x]=depth[p]+1;
	for (int y : adj[x]) if (y!=p) {
		if (y==1) z=x;  // back-edge from z to 1 (root)
		dfs(y,x);
	}
}

void go(int from, int to, VI &moves, VI &keys) {
	if (depth[from]<depth[to]) {
		VI movesR;
		go(to,from,movesR,keys);
		movesR.pop_back(); reverse(ALL(movesR)); movesR.push_back(to);
		moves.insert(moves.end(),ALL(movesR));
		return;
	}
	while (from!=to) {
		int x=prv[from];
		moves.push_back(x);
		keys.push_back(edge[{x,from}]);
		from=x;
	}
}

void print(VI kA, VI mA1, VI dA, VI mA2, VI kB, VI mB1, VI mB2) {
	for (int x : kA) cout << x << ' '; cout << endl;
	for (int x : kB) cout << x << ' '; cout << endl;
	for (int x : mA1) cout << "MOVE " << x << '\n';
	cout << "DROP"; for (int x : dA) cout << ' ' << x; cout << "\n";
	for (int x : mA2) cout << "MOVE " << x << '\n';
	cout << "DONE\n";
	for (int x : mB1) cout << "MOVE " << x << '\n';
	cout << "GRAB\n";
	for (int x : mB2) cout << "MOVE " << x << '\n';
	cout << "DONE\n";
}

int main() {
	//freopen("test.in","r",stdin);
	cin >> n >> m;
	REP (i,m) {
		int x,y;
		cin >> x >> y;
		adj[x].push_back(y); adj[y].push_back(x);
		edge[{x,y}]=i; edge[{y,x}]=i;
	}
	REP (i,n) reverse(ALL(adj[i]));
	dfs(1);
	if (z==-1) {
		cerr << "case 0" << endl;
		cout << "No solution" << endl;
		return 0;
	}
	assert(vis[0]);
	int x=z;
	for (az[x]=1;x!=1;x=prv[x]) az[prv[x]]=1;
	x=0;
	while (!az[x]) x=prv[x];
	VI emp;
	if (az[0]) {  // no drop, 0->1, 1-z->0
		cerr << "case 1" << endl;
		VI mA, kA;
		go(0,1,mA, kA);
		VI mB={z}, kB={edge[{1,z}]};
		go(z,0,mB, kB);
		print(kA,emp,emp,mA,kB,emp,mB);
	} else if (x==1) {  // 0->1->z(drop)->1, 1-z(grab)-1->0
		cerr << "case 2" << endl;
		VI mC, kC;
		go(0,1,mC,kC);
		VI mA1=mC, kA=kC, mA2;
		go(1,z,mA1,kA);
		go(z,1,mA2,emp);
		VI mB1={z}, kB={edge[{1,z}]}, mB2={1};
		go(1,0,mB2,emp);
		print(kA,mA1,kC,mA2,kB,mB1,mB2);
	} else {  // 0->x(drop)->1, 1-z->x(grab)->0
		cerr << "case 3" << endl;
		VI mC, kC;
		go(0,x,mC,kC);
		VI mA1=mC, kA=kC, mA2;
		go(x,1,mA2,kA);
		VI mB1={z}, kB={edge[{1,z}]}, mB2;
		go(z,x,mB1,kB);
		go(x,0,mB2,emp);
		print(kA,mA1,kC,mA2,kB,mB1,mB2);
	}
	return 0;
}
