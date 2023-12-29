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

#define maxn 100000

using namespace std;
typedef long long ll;

const ll MOD = 1000000007LL;
const float PI =  3.1415926535897932384626433;
double xa, ya, xb, yb, xc, yc, r;

double hyp(double x, double y){
    // math hypot has lower precision
    return sqrt(x*x+y*y);
}

double comp_dist(double ang){
    double x = xc+sin(ang)*r;
    double y = yc+cos(ang)*r;
    return hyp(xa-x, ya-y)+hyp(xb-x, yb-y);
}

void solve(){
    scanf("%lf%lf%lf%lf%lf%lf%lf",&xa,&ya,&xb,&yb,&xc,&yc,&r);
    if(hyp(xa-xc,ya-yc) <= r || hyp(xb-xc,yb-yc) <= r){
        printf("%.12lf\n",hyp(xa-xb,ya-yb));
        return;
    }
    double lang = 0, rang = 2*PI;
    FOR(cnt1,100){
        double best_ang=lang, best_d=1e9;
        int imax = cnt1==0?5000:10;
        double margin = (rang-lang)/imax;
        FOR(i,imax){
            double ang = lang+i*margin;
            double d = comp_dist(ang);
            if (d<best_d){best_ang=ang; best_d=d;}
        }
        lang = best_ang-margin;
        rang = best_ang+margin;
    }
    printf("%.12lf\n",comp_dist(lang));
}

int main(){
    int t;
    scanf("%d",&t);
    while(t--)solve();
	return 0;
}
