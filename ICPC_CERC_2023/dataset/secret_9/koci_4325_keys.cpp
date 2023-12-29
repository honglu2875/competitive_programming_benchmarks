#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cstring>

#define maxn 100010
using namespace std;

int n,m;

vector<int> neigh[maxn];
map<pair<int, int>, int> edge_ids;
int p_cyc[maxn];
int p_path[maxn];

int get_id(int x, int y){
    if (y<x)swap(x,y);
    return edge_ids[make_pair(x,y)];
}

void dfs_cyc(int node, int parent){
    // finds cycle from 1 to 1
    p_cyc[node]=parent;
    if(node==1 && parent != -1)return;
    for(int i=0;i<neigh[node].size();i++){
        int child = neigh[node][i];
        if (p_cyc[child]!=-1 || child==parent)continue;
        dfs_cyc(child,node);
    }
}

void dfs_path(int node){
    if(node==0)return;
    for(int i=0;i<neigh[node].size();i++){
        int child = neigh[node][i];
        if (p_path[child]!=-1)continue;
        p_path[child]=node;
        dfs_path(child);
    }
}

int main(){
    // first, a travels from 0 to 1,
    // then, b travels from 1 to 0
    scanf("%d%d", &n,&m);
    for(int i=0;i<m;i++){
        int a,b;
        scanf("%d%d",&a,&b);
        neigh[a].push_back(b);
        neigh[b].push_back(a);
        if(a>b)swap(a,b);
        edge_ids[make_pair(a,b)]=i;
    }
    memset(p_cyc, -1, sizeof(p_cyc));
    memset(p_path, -1, sizeof(p_path));
    dfs_cyc(1, -1);
    if(p_cyc[1]==-1){
        // 1 is not part of a cycle
        printf("No solution\n");
        return 0;
    }
    dfs_path(1);
    vector<int> nodes_in_path;
    vector<int> nodes_in_cycle;
    set<int> cyc_set;
    // make path a vector
    int node = 0;
    while(node != 1){
        nodes_in_path.push_back(node);
        node=p_path[node];
    }
    nodes_in_path.push_back(1);
    reverse(nodes_in_path.begin(), nodes_in_path.end());
    // make cycle a vector
    node = p_cyc[1];
    while(node != 1){
        nodes_in_cycle.push_back(node);
        cyc_set.insert(node);
        node=p_cyc[node];
    }
    nodes_in_cycle.push_back(node);
    reverse(nodes_in_cycle.begin(), nodes_in_cycle.end());
    // find the latest node in path that appears in the cycle
    int last_common_node=1;
    for(int i=0; i<nodes_in_path.size(); i++){
        int curr = nodes_in_path[i];
        if(cyc_set.find(curr)!=cyc_set.end())last_common_node=curr;
    }
    // case analysis:
    // Case 1: 0 lies on the cycle:
    //       A goes from 0 to 1 one way, B 1 to 0 following the other side of the cycle
    //       No key exchange required
    //
    // Case 2: path and cycle disjoint:
    //       A goes from 0 to 1, then the long way to p_cyc[1], drops keys for the path, returns the long way to 1.
    //       B goes from 1 to p_cyc[1] (direct), picks up keys, returns to 1, follows path 1->0
    //
    // Case 3: path and cycle overlap:
    //       Find the first node in 0->1 that exists in the cycle, let's call it X
    //       A goes from 0 to X, drops the used keys, then follows one of the cycle routes to 1.
    //       B takes the other route to 1->X, picks up the keys, then goes X to 0
    vector<int> shared_keys;
    for(int i=nodes_in_path.size()-1; nodes_in_path[i]!=last_common_node; i--){
        shared_keys.push_back(get_id(nodes_in_path[i],nodes_in_path[i-1]));
    }
    if (last_common_node == 0){
        vector<int> a_keys;
        vector<int> b_keys;
        int zero_idx = -1;
        for(int i=0;i<nodes_in_cycle.size()-1;i++){
            int curr_key = get_id(nodes_in_cycle[i], nodes_in_cycle[i+1]);
            if (nodes_in_cycle[i]==0)zero_idx=i;
            if(zero_idx == -1)a_keys.push_back(curr_key);
            else b_keys.push_back(curr_key);
        }
        b_keys.push_back(get_id(nodes_in_cycle.back(),1));
        for(int i=0;i<a_keys.size();i++)
            printf("%d%c",a_keys[i],i==a_keys.size()-1?'\n':' ');
        for(int i=0;i<b_keys.size();i++)
            printf("%d%c",b_keys[i],i==b_keys.size()-1?'\n':' ');
        for(int i=zero_idx-1;i>=0;i--)
            printf("MOVE %d\n",nodes_in_cycle[i]);
        printf("DONE\n");
        for(int i=nodes_in_cycle.size()-1;i>=zero_idx;i--)
            printf("MOVE %d\n",nodes_in_cycle[i]);
        printf("DONE\n");
    } else if (last_common_node == 1){
        vector<int> a_keys;
        for(int i=0; i<nodes_in_cycle.size()-1;i++)
            a_keys.push_back(get_id(nodes_in_cycle[i],nodes_in_cycle[i+1]));
        // A has starting keys shared_keys and a_keys
        for(int i=0;i<shared_keys.size();i++){
            if(i!=0)printf(" ");
            printf("%d",shared_keys[i]);
        }
        for(int i=0;i<a_keys.size();i++)
            printf(" %d",a_keys[i]);
        printf("\n");
        // B starts with key between 1 and the last node in the cycle
        printf("%d\n",get_id(1,nodes_in_cycle.back()));
        for(int i=nodes_in_path.size()-2;i>=0;i--)
            printf("MOVE %d\n",nodes_in_path[i]);
        for(int i=1; i<nodes_in_cycle.size();i++)
            printf("MOVE %d\n",nodes_in_cycle[i]);
        printf("DROP");
        for(int i=0;i<shared_keys.size();i++)printf(" %d",shared_keys[i]);
        printf("\n");
        for(int i=nodes_in_cycle.size()-2;i>=0;i--)
            printf("MOVE %d\n",nodes_in_cycle[i]);
        printf("DONE\n");
        printf("MOVE %d\n",nodes_in_cycle.back());
        printf("GRAB\n");
        for(int i=0;i<nodes_in_path.size();i++)
            printf("MOVE %d\n",nodes_in_path[i]);
        printf("DONE\n");
    } else {
        int x_idx;
        vector<int> a_keys;
        vector<int> b_keys;
        for(x_idx=0;nodes_in_cycle[x_idx]!=last_common_node;x_idx++)
            a_keys.push_back(get_id(nodes_in_cycle[x_idx],nodes_in_cycle[x_idx+1]));
        for(int i=x_idx;i<nodes_in_cycle.size()-1;i++)
            b_keys.push_back(get_id(nodes_in_cycle[i],nodes_in_cycle[i+1]));
        b_keys.push_back(get_id(nodes_in_cycle.back(),1));
        // A has starting keys shared_keys and a_keys
        for(int i=0;i<shared_keys.size();i++){
            if(i!=0)printf(" ");
            printf("%d",shared_keys[i]);
        }
        for(int i=0;i<a_keys.size();i++)
            printf(" %d",a_keys[i]);
        printf("\n");
        for(int i=0;i<b_keys.size();i++){
            if(i!=0)printf(" ");
            printf("%d",b_keys[i]);
        }
        printf("\n");
        // A moves to the cycle
        for(int i=nodes_in_path.size()-2;nodes_in_path[i]!=last_common_node;i--)
            printf("MOVE %d\n",nodes_in_path[i]);
        printf("MOVE %d\n",last_common_node);
        printf("DROP");
        for(int i=0;i<shared_keys.size();i++)printf(" %d",shared_keys[i]);
        printf("\n");
        for(int i=x_idx-1;i>=0;i--)
            printf("MOVE %d\n",nodes_in_cycle[i]);
        printf("DONE\n");
        for(int i=nodes_in_cycle.size()-1;i>=x_idx;i--)
            printf("MOVE %d\n",nodes_in_cycle[i]);
        printf("GRAB\n");
        int i;
        for(i=0;nodes_in_path[i]!=last_common_node;i++);
        i++;
        for(;i<nodes_in_path.size();i++)printf("MOVE %d\n",nodes_in_path[i]);
        printf("DONE\n");
    }
    return 0;
}
