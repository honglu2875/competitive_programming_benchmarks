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
template<class A, class B> ostream& operator<<(ostream& os, pair<A,B> &p) {
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

#define INF 1000000001
#define N 300000
#define K 500

struct Interval {
	int s,e;  // start, end
	int b;  // bucket
	int a;  // active
	int bi;  // index within a bucket
};

struct Bucket {
	int n;  // size
	int ints[K];  // bucket interval ids sorted by starts
	int starts[K];  // list of start value
	PII min_end[K+1];  // (end, ind)
	int fin[K+1];  // finish index
	int cnt[K+1];  // number of moves
};

int ops[N+1], T;

Interval ints[N+1];
PII rights[N+1];
int n;

int nb;
Bucket buckets[N/K+1];

int pos[N+1];  // interval with the last cut
int cnt[N+1];  // number of cuts

int gt[N+1];  // first > in a bucket

int main() {
	//freopen("test.in","r",stdin);
	//freopen("test.out","w",stdout);
	// read operations
	int x,y;
	scanf("%d",&T);
	ints[0]={-1,-1}; rights[0]={-1,0};
	FOR (t,1,T) {
		scanf("%d",&x);
		if (x<0) ops[t]=x;
		else {
			scanf("%d",&y);
			n++;
			ints[n]={x,y};
			rights[n]={y,n};
			ops[t] = n;
		}
	}
	// organize into buckets
	sort(rights,rights+n+1);
	nb = n/K+1;
	FOR (i,0,n) {
		auto [e,ind] = rights[i];
		int bi=i/K;
		ints[ind].b = bi;
		auto &b = buckets[bi];
		b.starts[b.n] = ints[ind].s;
		b.ints[b.n] = ind;
		b.n++;
	}
	REP (bi,nb) {
		auto &b = buckets[bi];
		sort(b.starts, b.starts+b.n);
		sort(b.ints, b.ints+b.n, [](int i, int j) { return ints[i].s < ints[j].s; });
		REP (i,b.n) {
			int id=b.ints[i];
			ints[id].bi=i;
		}
	}
	// process buckets
	FOR (t,1,T) { pos[t]=0; cnt[t]=0; }
	REP (bi,nb) {
		auto &b = buckets[bi];
		// init
		REP (i,b.n+1) { b.min_end[i]={INF,-1}; b.fin[i]=i; b.cnt[i]=0; }
		// compute locations within current bucket for all end points
		int j=0;
		FOR (ei,0,n) {
			auto [x,i]=rights[ei];
			while (j<b.n && b.starts[j]<=x) j++;
			gt[i]=j;
		}
		// process all times
		FOR (t,1,T) {
			// (de)activate if necessary (change in current bucket)
			int ind=abs(ops[t]);
			if (ints[ind].b==bi) {
				ints[ind].a=ops[t]>0;
				// update bucket
				for (int i=b.n-1;i>=0;i--) {
					// update minimum ends
					int id=b.ints[i];
					if (!ints[id].a) b.min_end[i]=b.min_end[i+1];
					else b.min_end[i]=min(b.min_end[i+1], make_pair(ints[id].e, id));
					// update all finish positions and counts
					b.fin[i]=i; b.cnt[i]=0;
					int j=gt[id];
					int id2=b.min_end[j].second;
					if (id2!=-1) {
						int i2=ints[id2].bi;
						b.fin[i]=b.fin[i2];
						b.cnt[i]=b.cnt[i2]+1;
					}
				}
			}
			// jump into bucket
			int id=pos[t];
			int j=gt[id];  // where does it land within the bucket
			int id2=b.min_end[j].second;
			if (id2!=-1) {  // something active exists in the bucket
				pos[t]=id2; cnt[t]++;
				int i2=ints[id2].bi;
				// jump within bucket
				cnt[t]+=b.cnt[i2];
				pos[t]=b.ints[b.fin[i2]];
			}
		}
	}
	FOR (t,1,T) printf("%d\n",cnt[t]);
	return 0;
}
