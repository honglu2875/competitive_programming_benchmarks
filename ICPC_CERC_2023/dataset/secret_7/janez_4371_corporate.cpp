#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cassert>
#include <random>
#include <algorithm>
#include <cstring>
#include <stack>
#include <utility>
#include <ctime>
#include <unordered_map>
using namespace std;

#define Assert assert
//#define Assert(x) 

struct Node
{
    string name;
    int parent = -1, lineNo = -1;
    vector<int> children;
    Node(string name_) : name(move(name_)) { }
};

vector<Node> nodes;
unordered_map<string, int> nameToIdx;

int NameToIdx(const string &s) 
{
    //fprintf(stderr, "s = \"%s\"\n", s.c_str());
    Assert(s.length() > 0); Assert(s.length() <= 10);
    for (char c : s) Assert(('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z')); 
    auto [it, isNew] = nameToIdx.try_emplace(s, (int) nodes.size());
    if (isNew) nodes.emplace_back(s);
    return it->second;
}

enum { MinEmployees = 2, MaxEmployees = 600, MaxB = 2048 };
enum { M = 1'000'000'007, ChkBits = 30 };

void Traverse(int curNode, vector<string> &names, string &bits)
{
    auto &N = nodes[curNode];
    names.emplace_back(N.name);
    for (int child : N.children)
    {
        bits += '0'; 
        Traverse(child, names, bits);
        bits += '1';
    }
}

unsigned int Checksum(const vector<string> &names, const string &bits)
{
    string s; for (const string name : names) { s += name; s += '#'; }
    s += bits; return (unsigned int) (hash<string>{}(s) % M);
}

template<typename R>
void PrepareOutputOrder(int curNode, R &r, vector<int> &dest)
{
    dest.emplace_back(curNode);
    vector<int> children = nodes[curNode].children;
    shuffle(children.begin(), children.end(), r);
    for (int child : children) PrepareOutputOrder(child, r, dest);
}

int main()
{
    //ifstream ifs("corporate\\binary_3.in"); cin.rdbuf(ifs.rdbuf());
    //ifstream ifs("corporate99b.in"); cin.rdbuf(ifs.rdbuf());
    string line; getline(cin, line);
    if (line == "ENCODE")
    {
        for (int lineNo = 0; ; ++lineNo)
        {
            line.clear(); getline(cin, line);
            if (line.empty()) break;
            int n = line.length(), i = 0;
            while (i < n && line[i] != ':') ++i;
            Assert(i < n);
            string parentName = line.substr(0, i); ++i;
            //printf("parent \"%s\"\n", parentName.c_str());
            int parentIdx = NameToIdx(parentName);
            Assert(nodes[parentIdx].lineNo < 0); nodes[parentIdx].lineNo = lineNo; // each parent must be listed only once
            Assert(i < n && line[i] == ' '); // each parent must have at least one child
            while (i < n)
            {
                Assert(line[i] == ' '); 
                while (i < n && line[i] == ' ') ++i;
                Assert(i < n);
                int j = i; while (j < n && line[j] != ' ') ++j;
                Assert(j > i); 
                string childName = line.substr(i, j - i); i = j;
                //printf("child \"%s\"\n", childName.c_str());
                int childIdx = NameToIdx(childName);
                nodes[parentIdx].children.emplace_back(childIdx);
                nodes[childIdx].parent = parentIdx;
            }
        }
        Assert(cin.eof());
        int n = nodes.size(); Assert(MinEmployees <= n); Assert(n <= MaxEmployees);
        const int root = 0;
        for (int i = 0; i < n; ++i)
        {
            // Parents must be listed before their children.
            int parent = nodes[i].parent;
            if (i == root) { Assert(parent < 0); continue; }
            // Parents must be listed before their children.
            Assert(nodes[parent].lineNo >= 0); 
            if (nodes[i].lineNo >= 0) Assert(nodes[i].lineNo > nodes[parent].lineNo);
        }
        vector<string> nameList; string bits;
        Traverse(root, nameList, bits);
        Assert(bits.length() == (n - 1) * 2);
        // Add some random bits and a checksum.
        mt19937_64 r(123); //random_device r; 
        while (bits.length() + ChkBits < MaxB) bits += '0' + uniform_int_distribution<int>(0, 1)(r);
        unsigned int checksum = Checksum(nameList, bits);
        for (int b = 0; b < ChkBits; ++b) bits += '0' + ((checksum >> b) & 1);
        mt19937_64 r2(checksum); for (int i = 1; i < nameList.size(); ++i) swap(nameList[i], nameList[uniform_int_distribution<int>(0, i)(r2)]);
        for (const string &name : nameList) cout << name << endl;
        cout << bits << endl;
    }
    else if (line == "DECODE")
    {
        // Read the data.
        vector<string> nameList;
        while (true)
        {
            string s; getline(cin, s); 
            if (s.empty()) { Assert(cin.eof()); break; }
            nameList.emplace_back(s);
        }
        Assert(nameList.size() >= MinEmployees + 1);
        string bits = nameList.back(); nameList.pop_back();
        const int n = nameList.size(); Assert(n <= MaxEmployees);
        Assert(bits.length() == MaxB); for (char c : bits) Assert(c == '0' || c == '1');
        // Check the checksum.
        unsigned int checksum = 0; 
        for (int i = 0; i < ChkBits; ++i) checksum |= ((unsigned int) (bits[bits.length() - ChkBits + i] - '0')) << i;
        mt19937_64 r2(checksum); vector<int> swapList;
        for (int i = 1; i < nameList.size(); ++i)
            swapList.emplace_back(uniform_int_distribution<int>(0, i)(r2));
        for (int i = 1; i < nameList.size(); ++i) 
            swap(nameList[nameList.size() - i], nameList[swapList[nameList.size() - i - 1]]);
        unsigned int checksum2 = Checksum(nameList, bits.substr(0, bits.length() - ChkBits));
        Assert(checksum == checksum2);
        // Reconstruct the tree.
        int root = NameToIdx(nameList[0]); Assert(root == 0);
        int curNode = root, nextName = 1;
        for (int i = 0; i < 2 * (n - 1); ++i)
        {
            if (bits[i] == '0') // move down
            {
                Assert(nextName < nameList.size());
                int childNode = NameToIdx(nameList[nextName++]);
                Assert(childNode == nodes.size() - 1);
                nodes[curNode].children.emplace_back(childNode);
                nodes[childNode].parent = curNode;
                curNode = childNode;
            }
            else // move up
            {
                Assert(bits[i] == '1');
                curNode = nodes[curNode].parent; Assert(curNode >= 0);
            }
        }
        Assert(curNode == root); Assert(nextName == nameList.size());
        mt19937_64 r(123); //random_device r; 
        vector<int> outputOrder; PrepareOutputOrder(root, r, outputOrder);
        Assert(outputOrder.size() == n);
        for (int lineNo = 0; lineNo < outputOrder.size(); ++lineNo)
        {
            int nodeIdx = outputOrder[lineNo]; auto &N = nodes[nodeIdx];
            Assert(N.lineNo < 0); N.lineNo = lineNo;
            if (N.children.empty()) continue;
            cout << N.name << ':';
            for (int child : N.children) cout << ' ' << nodes[child].name;
            cout << endl;
        }
    }
    else Assert(false);
    return 0;
}