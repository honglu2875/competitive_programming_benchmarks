#include <iostream>
#include <iomanip>
#include <climits>
#include <stack>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <deque>
#include <queue>
#include <set>
#include <map>
#include <cassert>

#define FOR(i,n) for(int i=0,_n=n;i<_n;i++)
#define FORR(i,s,n) for(int i=s,_n=n;i<_n;i++)
#define mp make_pair
#define pb push_back
#define pii pair<int,int>
#define pli pair<ll,int>
#define vi vector<int>
#define fs first
#define sec second

#define maxn 1000

using namespace std;
typedef long long ll;

const ll INF=4e18L;

struct edge{
	int a,b;
    ll c,flow;
};

int n,s,t;//number, source, sink
int d[maxn];
int ptr[maxn];
int q[maxn];

vector <edge> e;
vector <int> graf[maxn];

void add_edge(int a, int b, ll c){
	edge e1 = {a,b,c,0};
	edge e2 = {b,a,0,0};
	graf[a].pb((int)e.size());
	e.pb(e1);
	graf[b].pb((int)e.size());
	e.pb(e2);
}

bool bfs(){
	int qh=0,qt=0;
	q[qt++]=s;
	memset(d,-1,n*sizeof(d[0]));
	d[s]=0;
	while(qh<qt&&d[t]==-1){
		int tr=q[qh++];
		FOR(i,graf[tr].size()){
			int id=graf[tr][i];
			int to=e[id].b;
			if(d[to]==-1 && e[id].flow<e[id].c){
				q[qt++]=to;
				d[to]=d[tr]+1;
			}
		}
	}
	return d[t]!=-1;
}

ll dfs(int v, ll flow){
	if(flow==0)return 0;
	if(v==t)return flow;
	for(;ptr[v]<(int)graf[v].size();ptr[v]++){
		int id=graf[v][ptr[v]];
		int to=e[id].b;
		if(d[to]!=d[v]+1)continue;
		ll pushed=dfs(to,min(flow,e[id].c-e[id].flow));
		if(pushed!=0){
			e[id].flow+=pushed;
			e[id^1].flow-=pushed;
			return pushed;
		}
	}
	return 0;
}

ll dinic(){
	ll flow=0;
	while(1){
		if(!bfs())break;
		memset(ptr,0,n*sizeof(ptr[0]));
		while(ll pushed=dfs(s,INF))flow+=pushed;
	}
	return flow;
}

ll ingredients[maxn];
ll tools[maxn];
ll cake_prices[maxn];

int main(){
    int ni,nc,nt;
    scanf("%d%d%d",&ni,&nc,&nt);
    n=2+nc+nt, s=0, t=1+nc+nt;
    FOR(i,nc)scanf("%lld",cake_prices+i);
    FOR(i,ni)scanf("%lld",ingredients+i);
    FOR(i,nc)assert(ingredients[i]<=1e8);
    FOR(i,nt){
        scanf("%lld",tools+i);
        add_edge(nc+1+i,t,tools[i]);
    }
    ll total_profit=0;
    FOR(i,nc){
        ll cost=0;
        FOR(j,ni){
            ll amount;
            scanf("%lld",&amount);
            assert(amount<=1e8);
            cost+=amount*ingredients[j];
        }
        if (cost>cake_prices[i])continue;
        add_edge(s,1+i, cake_prices[i]-cost);
        total_profit+=cake_prices[i]-cost;
    }
    FOR(i,nc){
        int n_tools=0;
        scanf("%d",&n_tools);
        FOR(j,n_tools){
            int tool_id;
            scanf("%d",&tool_id); //1-based indexing
            add_edge(i+1,tool_id+nc, INF);
        }
    }
    printf("%lld\n",total_profit-dinic());
	return 0;
}
