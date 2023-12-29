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

int main() {
	string line="";
	map<string,int> f;
	vector<pair<int,string>> todo = {{-1,"------"}, {1,"======"}};
	for (auto [d, fin] : todo) {
		while (1) {
			getline(cin,line);
			if (line==fin) break;
			stringstream ss(line);
			int s,e;
			string u;
			ss >> s >> e >> u;
			f[u]+=d*(e-s);
		}
	}
	int cnt=0;
	for (auto [user, dif] : f) if (dif!=0) {
		cout << user << " ";
		if (dif>0) cout << "+" << dif << endl;
		else cout << dif << endl;
		cnt++;
	}
	if (cnt==0) cout << "No differences found." << endl;
	return 0;
}
