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
#include <map>
using namespace std;

//#define Assert assert
#define Assert(x) 

map<string, int> differences;

enum { MaxTime = 1000, MinNameLen = 3, MaxNameLen = 20 };

void ReadSchedule(int coef, const string &endLine)
{
    string line;
    int nLines = 0, curTime = 0;
    while (true)
    {
        line.clear(); getline(cin, line);
        if (line == endLine) break;
        int i = 0, n = line.length();
        int startTime = -1;
        while (i < n && line[i] >= '0' && line[i] <= '9')
        {
            startTime = (startTime < 0 ? 0 : 10 * startTime) + (line[i++] - '0'); 
            Assert(startTime <= MaxTime);
        }
        Assert(startTime == curTime); 
        Assert(i < n); Assert(line[i] == ' '); ++i;
        int endTime = -1;
        while (i < n && line[i] >= '0' && line[i] <= '9')
        {
            endTime = (endTime < 0 ? 0 : 10 * endTime) + (line[i++] - '0'); 
            Assert(endTime <= MaxTime);
        }
        Assert(endTime > startTime); 
        Assert(i < n); Assert(line[i] == ' '); ++i;
        string name = line.substr(i);
        for (char c : name) Assert(c >= 'a' && c <= 'z');
        Assert(name.length() >= MinNameLen); Assert(name.length() <= MinNameLen);
        curTime = endTime;
        auto [it, isNew] = differences.try_emplace(name, 0);
        it->second += (endTime - startTime) * coef; ++nLines;
        //printf(" - %s %d..%d -> now %d\n", name.c_str(), startTime, endTime, it->second);
    }
    Assert(nLines > 0);
}

int main()
{
    ReadSchedule(-1, "------");
    ReadSchedule(1, "======");
    string line; getline(cin, line); Assert(line.empty()); Assert(cin.eof());
    bool any = false;
    for (const auto& [name, difference] : differences) if (difference != 0) printf("%s %+d\n", name.c_str(), difference), any = true;
    if (! any) printf("No differences found.\n");
    return 0;
}