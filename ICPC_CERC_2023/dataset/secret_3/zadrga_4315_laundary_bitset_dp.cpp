#include <bits/stdc++.h>
 
using namespace std;

#define maxn 30011
#define maxl 300011
#define INF 1000 * 1000 * 1000 + 9

typedef long long ll;
typedef long double ld;
typedef pair<int, int> pii;

class Sheet{
    public:
        int d, t_short, t_long;

    Sheet(int d, int t_short, int t_long){
        this->d = d;
        this->t_short = t_short;
        this->t_long = t_long;
    }

    bool operator< (const Sheet& s){
        return this->t_long < s.t_long;
    }
};

/*
dp[i][j]: whether it is possible to hang the first i items and use j meters on the first line
          (the space left on the second line can be directly deduced)
*/
bitset<maxl> dp[2];
bitset<maxl> dp_reverse[2];
int ans[maxl];
vector<Sheet> v;
vector<int> total_length;


int main(){
    int N, Q;
    scanf("%d%d", &N, &Q);

    assert(1 <= N && N <= 30000);
    assert(1 <= Q && Q <= 300000);

    v.push_back(Sheet(0, 0, 0));
    int max_t_short = -1;
    for(int i = 0; i < N; i++){
        int d, t_short, t_long;
        scanf("%d%d%d", &d, &t_short, &t_long);
        v.push_back(Sheet(d, t_short, t_long));

        max_t_short = max(max_t_short, t_short);
        
        assert (1 <= d && d <= 300000);
        assert (1 <= t_short && t_short <= 1000*1000*1000);
        assert (1 <= t_long && t_long <= 1000*1000*1000);
        assert (t_short <= t_long);
    }

    sort(v.begin(), v.end());

    for(int i = 0; i < v.size(); i++){
        if (i == 0){
            total_length.push_back(v[i].d);
        }
        else{
            total_length.push_back(total_length[i - 1] + v[i].d);
        }
    }

    int last_L = maxl - 1;

    // All sheets can be hung over both lines.
    while(last_L >= 0 && last_L >= total_length[N]){
        ans[last_L] = max_t_short;
        last_L--;
    }


    dp[0][maxl - 1] = 1;
    dp_reverse[0][0] = 1;
    for(int i = 1; i <= N; i++){
        int cur = i % 2;
        int last = 1 - cur;


        dp[cur] = dp[last] | (dp[last] >> v[i].d);
        dp_reverse[cur] = dp_reverse[last] | (dp_reverse[last] << v[i].d);

        // Let i be the first sheet that is hung over only a single line,
        // while sheets (i + 1), ..., N are hung over both lines. 

        // If it is possible to dry all the clothes, at least maximum of 
        // all short drying times is required.

        int T = max(max_t_short, v[i].t_long);
        int both_lines = total_length[N] - total_length[i];

        // Determine the shortest required length of lines so that the 
        // sheets 1, ..., i can be hung over only a single line.
        // The optimal value is one of the two closest to total_length[i] / 2.
        
        int min_L_single = INF;
        int mid = total_length[i] / 2 - 1;

        int next_L_single = (maxl - 1) - dp[cur]._Find_next((maxl - 1) - mid);
        min_L_single = min(min_L_single, max(next_L_single, total_length[i] - next_L_single));


        int prev_L_single = dp_reverse[cur]._Find_next(mid);
        min_L_single = min(min_L_single, max(prev_L_single, total_length[i] - prev_L_single));


        int min_L_required = min_L_single + both_lines;
        
        // Update answers. 
        while(last_L >= 0 && last_L >= min_L_required){
            ans[last_L] = T;
            last_L--;
        }     
    }

    // It is impossible to hang the sheets.
    while(last_L >= 0){
        ans[last_L] = -1;
        last_L--;
    }

    for(int i = 0; i < Q; i++){
        int L;
        scanf("%d", &L);

        assert (1 <= L && L <= 300000);
        
        printf("%d\n", ans[L]);
    }

    return 0;
}