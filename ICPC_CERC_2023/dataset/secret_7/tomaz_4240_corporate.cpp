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

string mode;
map<string,VS> f;

string encoding="";
void encode(string name) {
	cout << name << endl;
	encoding+="0";
	for (string child : f[name]) encode(child);
	encoding+="1";
}

vector<string> names;
string line;
int line_ind=0, name_ind=0;
void decode() {
	string name=names[name_ind];
	name_ind++; line_ind++;
	while (line[line_ind]=='0') {
		f[name].push_back(names[name_ind]);
		decode();
	}
	line_ind++;
}

void print(string name) {
	if (!f[name].empty()) {
		cout << name << ":";
		for (string child : f[name]) cout << " " << child;
		cout << endl;
	}
	for (string child : f[name]) print(child);
}

int main() {
	//freopen("test_encode.in","r",stdin);
	//freopen("test_decode.in","r",stdin);
	getline(cin,mode);
	if (mode=="ENCODE") {
		string root="";
		while (getline(cin,line)) {
			stringstream ss(line);
			string parent,child;
			ss >> parent;
			parent.pop_back();
			if (root.empty()) root=parent;
			while (ss >> child) f[parent].push_back(child);
		}
		encode(root);
		cout << encoding << endl;
	}
	if (mode=="DECODE") {
		string root="", prev="";
		while (getline(cin,line)) {
			if (root.empty()) root=line;
			if (!prev.empty()) names.push_back(prev);
			prev=line;
		}
		decode();
		print(root);
	}
	return 0;
}
