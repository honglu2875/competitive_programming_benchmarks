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
#define DEBUG(x) { cerr << #x << " = "; cerr << x << endl; }
#define DEBUG_ITER(x) { cerr << #x << " = ["; for (auto _ : x) cerr << _ << "| "; cerr << "]" << endl; }
template<class A, class B> ostream& operator<<(ostream& os, pair<A,B> p) {
	os << "(" << p.first << ", " << p.second << ")";
	return os;
}

typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int,int> PII;
typedef vector<int> VI;
typedef vector<PII> VII;
typedef array<int,3> III;
typedef vector<string> VS;

#define N 100000
#define K 100000
#define MOD 1000000007

int n,m,k;
vector<int> adj0[N], adj[N], adj2[N];
set<int> adj_set[N];

int remaining, act[N];
int root;

int other(VI sez, VI invalid) {
	for (auto x : sez) if (!IN(x,invalid)) return x;
	assert(false);
}

void delEdge(int a, int b) {
	adj[a].erase(find(ALL(adj[a]),b));
	adj[b].erase(find(ALL(adj[b]),a));
}
void addEdge(int a, int b) {
	adj[a].PB(b);
	adj[b].PB(a);
}
void delNode(int a) {
	while (!adj[a].empty()) delEdge(a,adj[a][0]);
	act[a]=0; remaining--;
}

int inQ[N];
queue<int> q3;
void enqueue(int x) {
	if (adj[x].size()==3 && !inQ[x]) { q3.push(x); inQ[x]=1; }
}
int dequeue() {
	int x=q3.front(); q3.pop(); inQ[x]=0;
	return x;
}

int out[N];
int decompose() {
	if (remaining==4) {
		out[root]=0;
		REP (x,n) if (act[x] && x!=root) { out[x]=1; }
		return 1;
	}
	while (!q3.empty()) {
		int x=dequeue();
		if (x==root || !act[x] || adj[x].size()!=3) continue;
		int a=adj[x][0], b=adj[x][1], c=adj[x][2];
		for (auto [y,z] : vector<PII>{{a,b},{a,c},{b,c}}) {  // x-y-z
			if (y==root || z==root || adj[y].size()!=3 || adj[z].size()!=3) continue;
			if (IN(z,adj[y])) {  // triangle
				int x1=other(adj[x],{y,z}), y1=other(adj[y],{x,z}), z1=other(adj[z],{x,y});
				if (x1==y1 || x1==z1 || y1==z1) continue;  // common neighbours
				// compress into x
				delNode(y); delNode(z); addEdge(x,y1); addEdge(x,z1);
				enqueue(x); enqueue(x1); enqueue(y1); enqueue(z1);
				int ret=decompose();
				remaining+=2;
				act[y]=1; act[z]=1;
				out[x]=out[x1]; out[y]=out[y1]; out[z]=out[z1];
				return ret;
			} else {  // path center
				int v=other(adj[x],{y,z});  // inner parent
				if (!IN(v,adj[y]) || !IN(v,adj[z])) continue;  // not common
				delNode(x); addEdge(y,z);
				enqueue(y); enqueue(z); enqueue(v);
				int ret=decompose();
				remaining+=1;
				act[x]=1;
				out[x]=1;
				return ret;
			}
		}
	}
	return 0;
}

int isLeaf[N];

int par[N];
void parents(int x, int px=-1) {
	par[x]=px;
	if (px!=-1) adj[x].erase(find(ALL(adj[x]),px));
	if (isLeaf[x]) return;
	for (int y : adj[x]) if (y!=px) {
		parents(y,x);
	}
}

int vis[N];
vector<int> cyc;

int climb(int x) {
	while (x!=root && !vis[x]) { vis[x]=1; x=par[x]; }
	return x;
}

vector<int> ch[N];

void add(int &r, int d) { r=(r+d)%MOD; }

PII norm(int c, int l, int r) {
	if (l==c) return {0,r!=c};
	else return {1,2*(r!=c)-(r==l)};
}

int f[N][2][3];
int g[N][3], gr[N];

void solve(int x) {
	if (isLeaf[x]) {
		f[x][0][0]=1;
	} else {
		for (int y : ch[x]) solve(y);
		REP (l,2) {
			int cntC,cntU,cntUV;
			int y=ch[x][0];
			// first subtree
			REP (r,3) {
				g[0][r]=0;
				REP (c,min(k,4)) if (c!=0) {
					auto [l1,r1] = norm(c,l,r);
					cntC=f[y][l1][r1];
					add(g[0][r], cntC);
				}
				if (k>4) add(g[0][r], ((int64)(k-4)*cntC)%MOD);
			}
			// other subtrees
			for (int i=1;i<(int)ch[x].size();i++) {
				int y=ch[x][i];
				REP (r,3) {
					g[i][r]=0;
					REP (c,min(k,4)) if (c!=0) {
						cntC=0;
						REP (u,min(k,5)) {
							cntU=0;
							auto [l1,u1] = norm(0,l,u);
							REP (v,min(k,6)) if (u!=v) {
								auto [v1,r1] = norm(c,v,r);
								cntUV=((int64)g[i-1][u1]*f[y][v1][r1])%MOD;
								add(cntU, cntUV);
							}
							if (k>6) { add(cntU,((int64)(k-6)*cntUV)%MOD); }
							add(cntC, cntU);
						}
						if (k>5) { add(cntC,((int64)(k-5)*cntU)%MOD); }
						add(g[i][r], cntC);
					}
					if (k>4) add(g[i][r],((int64)(k-4)*cntC)%MOD);
				}
			}
			// final result
			REP (r,3) {
				f[x][l][r]=g[ch[x].size()-1][r];
			}
		}
	}
}

int main() {
	//freopen("test.in","r",stdin);
	//freopen("test.out","w",stdout);
	scanf("%d %d %d\n",&n,&m,&k);
	assert(4<=n && n<=N); assert(m<=2*(n-1)); assert(1<=k && k<=K);
	REP (i,m) {
		int a,b;
		scanf("%d %d",&a,&b); a--; b--;
		adj[a].PB(b); adj[b].PB(a);
		adj_set[a].insert(b); adj_set[b].insert(a);
	}
	REP (i,n) act[i]=1;
	REP (i,n) adj0[i] = adj[i];
	// remove deg-2 nodes
	remaining=n;
	vector<int> deg2;
	REP (x,n) if (adj_set[x].size()==2) {
		deg2.push_back(x);
		auto it=adj_set[x].begin();
		int a=*it, b=*(++it);
		remaining--;
		adj_set[a].erase(x); adj_set[b].erase(x);
		adj_set[a].insert(b); adj_set[b].insert(a);
		act[x]=0; adj_set[x].clear();
	}
	REP (x,n) adj[x] = vector<int>(ALL(adj_set[x]));
	REP (x,n) adj2[x] = adj[x];
	// find candidate roots (inner nodes)
	root = -1;
	REP (x,n) if (act[x] && adj[x].size()>3) root=x;
	vector<int> roots;
	if (root!=-1) roots.PB(root);
	else {
		REP (x,n) if (act[x]) {
			roots.PB(x);
			assert(adj[x].size()==3);
			roots.insert(roots.end(), ALL(adj[x]));
			break;
		}
	}
	// try to decompose Halin graph into inner and outter nodes
	int success=0;
	for (int cand : roots) {
		root = cand;
		// init
		REP (x,n) out[x]=-1;
		REP (x,n) adj[x] = adj2[x];
		while (!q3.empty()) dequeue();
		REP (x,n) if (act[x] && adj[x].size()==3) enqueue(x);
		// decompose
		success=decompose();
		if (!success) continue;
		REP (x,n) isLeaf[x]=(out[x]==1)?1:0;
		REP (x,n) adj[x] = adj0[x];  // restore initial graph
		REP (x,n) par[x] = -2;
		parents(root);  // try to reconstruct the tree
		REP (x,n) if (act[x] && par[x]==-2) { success=0; break; }  // active should be visited
		for (int x : deg2) if (par[x]==-2) { success=0; break; }  // deg2 should be part of the tree
		if (success) break;
	}
	assert(success);
	// find order of leaves on the cycle
	int start;
	REP (x,n) if (isLeaf[x]) { start=x; break; }
	int x=start, prev=-1;
	while (1) {
		cyc.PB(x);
		for (int y : adj[x]) {
			if (y!=prev && isLeaf[y]) { prev=x; x=y; break; }
		}
		if (x==start) break;
	}
	// find two adjacent with LCA=root (first, last)
	INIT(vis,0);
	climb(cyc.back());
	REP (i,cyc.size()) {
		if (climb(cyc[i])==root) { rotate(cyc.begin(),cyc.begin()+i,cyc.end()); break; }
	}
	// reorder children in the tree
	INIT(vis,0);
	for (int l : cyc) {
		int x=l;
		while (x!=root && !vis[x]) {
			vis[x]=1;
			int p = par[x];
			ch[p].PB(x);
			x = p;
		}
	}
	// count colorings
	solve(root);
	int st=0;
	add(st, ((int64)k*(k-1)%MOD)*f[root][0][1]%MOD);  // c=l, r
	add(st, ((int64)k*(k-1)%MOD)*f[root][1][0]%MOD);  // c=r, l
	add(st, ((int64)k*(k-1)*(k-2)%MOD)*f[root][1][2]%MOD);  // c, l, r
	printf("%d\n",st);
	return 0;
}
