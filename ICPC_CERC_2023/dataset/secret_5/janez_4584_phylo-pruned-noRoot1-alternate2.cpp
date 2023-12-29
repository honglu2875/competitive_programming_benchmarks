#define _CRT_SECURE_NO_WARNINGS 
#include <string>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cassert>
#include <random>
#include <algorithm>
#include <cstring>
#include <utility>
#include <ctime>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <list>
#include <thread>
using namespace std;

#define Assert assert
//#define Assert(x) 

using VertNo = int;
using EdgeNo = int;
typedef long long int llint;

struct Vert
{
    int deg = 0;
    bool exists = true;
    EdgeNo firstEdge = -1, lastEdge = -1;
    VertNo mergedInto = -1;
    int mergeFlag = -1;
    VertNo parent = -1;
};

struct Edge
{
    bool exists = true;
    VertNo u = -1, v = -1;
    EdgeNo uPrev = -1, uNext = -1, vPrev = -1, vNext = -1;

    VertNo Other(VertNo w) const { Assert(u == w || v == w); return (u == w) ? v : u; }
    EdgeNo Next(VertNo w) const { Assert(u == w || v == w); return (u == w) ? uNext : vNext; }
    EdgeNo Prev(VertNo w) const { Assert(u == w || v == w); return (u == w) ? uPrev : vPrev; }
};

struct MergeInfo
{
    vector<VertNo> inVerts;
    VertNo outVert = -1;
    vector<pair<EdgeNo, EdgeNo>> edges; // (inEdge, outEdge) pairs; for edges that connect two merged vertices, outEdge = -1
};

template<> struct std::hash<pair<VertNo, VertNo>> {
    inline size_t operator() (pair<VertNo, VertNo> p) const {
        return (size_t(p.first + 1) << 32) | size_t(p.second + 1); }
};

struct Graph
{
    vector<Vert> verts;
    vector<Edge> edges;
    int nVerts = 0, nEdges = 0; // only existent ones
    unordered_map<pair<VertNo, VertNo>, EdgeNo> edgeMap;

    VertNo AddVert() {
        VertNo v = verts.size(); verts.push_back(Vert{}); ++nVerts; return v; }

protected:

    void AddEdgeToLinkedLists(EdgeNo e)
    {
        auto &E = edges[e]; Assert(E.exists); 
        const VertNo u = E.u, v = E.v; Assert(u >= 0); Assert(v >= 0);
        auto &U = verts[u], &V = verts[v]; Assert(U.exists); Assert(V.exists);
        ++U.deg; ++V.deg;
        E.uPrev = U.lastEdge; 
        if (U.lastEdge < 0) U.firstEdge = e;
        else { 
            auto &EE = edges[U.lastEdge]; if (EE.u == u) { Assert(EE.uNext < 0); EE.uNext = e; }
            else { Assert(EE.v == u); Assert(EE.vNext < 0); EE.vNext = e; } }
        U.lastEdge = e;
        E.vPrev = V.lastEdge;
        if (V.lastEdge < 0) V.firstEdge = e;
        else {
            auto &EE = edges[V.lastEdge]; if (EE.u == v) { Assert(EE.uNext < 0); EE.uNext = e; }
            else { Assert(EE.v == v); Assert(EE.vNext < 0); EE.vNext = e; } }
        V.lastEdge = e;
    }

public:

    EdgeNo AddEdge(VertNo u, VertNo v)
    {
        Assert(u != v);
        auto &U = verts[u], &V = verts[v]; Assert(U.exists); Assert(V.exists);
        EdgeNo e = edges.size(); edges.push_back(Edge{}); auto &E = edges[e];
        { auto [it, isNew] = edgeMap.try_emplace({min(u, v), max(u, v)}, e); Assert(isNew); }
        E.u = u; E.v = v; 
        AddEdgeToLinkedLists(e); ++nEdges;
        return e;
    }
    EdgeNo GetEdge(VertNo u, VertNo v) const { 
        auto it = edgeMap.find({min(u, v), max(u, v)});
        return (it == edgeMap.end()) ? -1 : it->second; }
    EdgeNo GetEdgeIfExists(VertNo u, VertNo v) const { 
        EdgeNo e = GetEdge(u, v); if (e < 0) return e;
        auto &E = edges[e]; return (E.exists) ? e : -1; }
    bool EdgeExists(VertNo u, VertNo v) const { return GetEdgeIfExists(u, v) >= 0; }
    EdgeNo FindInNeighList(VertNo u, VertNo v) const { // tries to find 'v' in the list of u's neighbours
        auto &U = verts[u]; Assert(U.exists);
        for (EdgeNo e = U.firstEdge; e >= 0; ) {
            auto &E = edges[e]; Assert(E.exists);
            if (E.Other(u) == v) return e;
            e = E.Next(u); } 
        return -1; }
    void DelEdge(EdgeNo e) 
    {
        auto &E = edges[e]; Assert(E.exists);
        const VertNo u = E.u, v = E.v;
        auto &U = verts[u], &V = verts[v]; Assert(U.exists); Assert(V.exists);
        if (E.uPrev < 0) { Assert(U.firstEdge == e); U.firstEdge = E.uNext; }
        else { auto &EE = edges[E.uPrev]; if (EE.u == u) EE.uNext = E.uNext; else { Assert(EE.v == u); EE.vNext = E.uNext; } }
        if (E.uNext < 0) { Assert(U.lastEdge == e); U.lastEdge = E.uPrev; }
        else { auto &EE = edges[E.uNext]; if (EE.u == u) EE.uPrev = E.uPrev; else { Assert(EE.v == u); EE.vPrev = E.uPrev; } }
        if (E.vPrev < 0) { Assert(V.firstEdge == e); V.firstEdge = E.vNext; }
        else { auto &EE = edges[E.vPrev]; if (EE.u == v) EE.uNext = E.vNext; else { Assert(EE.v == v); EE.vNext = E.vNext; } }
        if (E.vNext < 0) { Assert(V.lastEdge == e); V.lastEdge = E.vPrev; }
        else { auto &EE = edges[E.vNext]; if (EE.u == v) EE.uPrev = E.vPrev; else { Assert(EE.v == v); EE.vPrev = E.vPrev; } }
        --U.deg; --V.deg; --nEdges;
        E.exists = false; E.uPrev = -1; E.uNext = -1; E.vPrev = -1; E.vNext = -1;
    }
    void DelVert(VertNo u) 
    {
        auto &U = verts[u]; Assert(U.exists);
        while (U.firstEdge >= 0) DelEdge(U.firstEdge);
        Assert(U.lastEdge < 0); Assert(U.deg == 0);
        U.exists = false; --nVerts;
    }
    void UndelVert(VertNo u)
    {
        auto &U = verts[u]; Assert(! U.exists);
        Assert(! U.exists); Assert(U.deg == 0); Assert(U.firstEdge < 0); Assert(U.lastEdge < 0);
        U.mergedInto = -1;
        U.exists = true; ++nVerts;
    }
    void UndelEdge(EdgeNo e)
    {
        auto &E = edges[e]; Assert(! E.exists);
        Assert(E.u >= 0); Assert(E.v >= 0); Assert(E.uPrev < 0); Assert(E.uNext < 0); Assert(E.vPrev < 0); Assert(E.vNext < 0);
        E.exists = true; AddEdgeToLinkedLists(e); ++nEdges;
    }
    
    int mergeFlag = 0;
    
    void Merge(MergeInfo &mi) // the caller should provide mi.inVerts; Merge() fills the rest while performing the merge
    {
        ++mergeFlag; 
        Assert(mi.inVerts.size() > 1);
        // Create the merged vertex.
        mi.outVert = AddVert();
        // Flag the vertices which participate in this marge.
        for (VertNo u : mi.inVerts) { auto &U = verts[u]; Assert(U.exists); U.mergeFlag = mergeFlag; }
        // Convert the edges and delete the vertices.
        mi.edges.clear();
        for (VertNo u : mi.inVerts) 
        {
            auto &U = verts[u];
            while (U.firstEdge >= 0)
            {
                EdgeNo eIn = U.firstEdge; auto &E = edges[eIn]; Assert(E.exists);
                Assert(E.u == u || E.v == u);
                VertNo v = (E.u == u) ? E.v : E.u; // the other end of this edge
                auto &V = verts[v]; Assert(V.exists); 
                EdgeNo eOut = -1;
                if (V.mergeFlag != mergeFlag) {
                    eOut = GetEdge(mi.outVert, v);
                    if (eOut >= 0) { Assert(edges[eOut].exists); eOut = -1; } // don't add a parallel edge
                    else eOut = AddEdge(mi.outVert, v); }
                mi.edges.emplace_back(eIn, eOut);
                DelEdge(eIn);
            }
            Assert(U.deg == 0); Assert(U.firstEdge < 0); Assert(U.lastEdge < 0);
            U.mergedInto = mi.outVert; 
        }
        for (VertNo u : mi.inVerts) DelVert(u);
    }

    void Unmerge(const MergeInfo &mi)
    {
        ++mergeFlag;
        auto &Z = verts[mi.outVert]; Assert(Z.exists);
        for (VertNo u : mi.inVerts) {
            auto &U = verts[u]; Assert(! U.exists); Assert(U.mergedInto == mi.outVert);
            UndelVert(u); }
        for (auto [eIn, eOut] : mi.edges) {
            auto &EI = edges[eIn]; Assert(! EI.exists); 
            if (eOut >= 0) {
                auto &EO = edges[eOut]; Assert(EO.exists); Assert(EO.u == mi.outVert || EO.v == mi.outVert); 
                DelEdge(eOut); }
            UndelEdge(eIn); }
        Assert(Z.deg == 0); Assert(Z.firstEdge < 0); Assert(Z.lastEdge < 0);
        DelVert(mi.outVert);
    }

private:

    static inline void Nop() { }

    bool RecoverTree_TryRoot(VertNo provisionalRoot, const VertNo tVerts[4],
        const vector<MergeInfo> &miThrees, const vector<MergeInfo> &miTwos)
    {
        #define XAssert(b) if (! (b)) return false; else Nop()
        Assert(provisionalRoot >= 0);
        for (int i = 0; i < 4; ++i) {
            VertNo u = tVerts[i]; auto &U = verts[u];
            U.parent = (u == provisionalRoot) ? -1 : provisionalRoot; }
        // Now undo the merges and propagate the parent information suitably.
        for (int mergeNo = int(miThrees.size()) - 1; mergeNo >= 0; --mergeNo)
        {
            const auto &mi = miThrees[mergeNo];
            VertNo z = mi.outVert; auto &Z = verts[z]; Assert(Z.exists);
            VertNo parent = Z.parent; Assert(parent >= 0); Assert(FindInNeighList(z, parent) >= 0);
            Unmerge(mi); 
            if (mi.inVerts.size() == 2) 
            {
                // The two source vertices will get the same parent as the destination vertex.
                for (VertNo u : mi.inVerts) {
                    verts[u].parent = parent;
                    XAssert(FindInNeighList(u, parent) >= 0); }
            }
            else 
            {
                // Otherwise this is a triangular merge.  The merged node must have its parent
                // as a neighbour, and exactly one of the unmerged nodes must also have this
                // parent as a neighbour.  For that unmerged node, the parent will be the same as for
                // the merged node, while that unmerged node in turn will become the parent of the
                // other two unmerged nodes.
                Assert(mi.inVerts.size() == 3);
                VertNo newChild = -1;
                for (auto [eIn, eOut] : mi.edges) {
                    if (eOut < 0) continue; // those were the edges that actually constituted the triangle
                    auto &EI = edges[eIn], EO = edges[eOut]; 
                    if (! ((EO.u == z && EO.v == parent) || (EO.u == parent && EO.v == z))) continue;
                    Assert(EI.u == parent || EI.v == parent);
                    VertNo u = EI.Other(parent);
                    newChild = u; break; }
                Assert(newChild >= 0);
                for (VertNo u : mi.inVerts) {
                    VertNo uParent = (u == newChild) ? parent : newChild;
                    auto &U = verts[u]; U.parent = uParent; if (uParent >= 0) Assert(FindInNeighList(u, uParent) >= 0); }
            }
            // If any of z's neighbours had z as a parent, update this to the appropriate
            // neighbour from the two or three nodes that had been merged to produce z.
            for (VertNo u : mi.inVerts) for (EdgeNo e = verts[u].firstEdge; e >= 0; ) { 
                auto &E = edges[e]; VertNo v = E.Other(u); auto &V = verts[v];
                if (V.parent == z) V.parent = u; 
                e = E.Next(u); }
        }
        // Lastly, undo the merges of 2-degree nodes.  
        for (int mergeNo = int(miTwos.size()) - 1; mergeNo >= 0; --mergeNo)
        {
            auto &mi = miTwos[mergeNo]; VertNo z = mi.outVert; auto &Z = verts[z];
            VertNo parent = Z.parent; Assert(mi.inVerts.size() == 2);
            VertNo u = mi.inVerts[0], v = mi.inVerts[1]; 
            if (parent >= 0)
            {
                // Remember that z was merged from two nodes, u and v, of which one had degree 2.
                // So the situation must be:
                //  After merge:                             Before merge:
                //       z ---------> parent                 u ----> v ------------> parent       OR    u ----> v --------------> parent
                //     ^^^^                                ^^^^                                         ^      ^^^
                //     ||||                                ||||                                         |      |||
                //     z's children*                    z's children*                              one of      z's other children
                //                                                                           z's children
                //     [*Actually those aren't necessarily just children.  z and u could be leaves
                //      with two connections therefore going to the previous and next leaf on the leaf cycle.]
                // Before the merge, all edges incident on z pointed into z, except for the one that pointed from z to its parent.
                // Now the direction of all those edges must remain the same, even though the endpoint z is replaced by one of u and v.
                // Moreover, the new edge between u and v must be directed towards the node which will then point to z's former parent.
                // We'll ensure that u is the one that points to z's former parent, so that v can then point to u.
                Assert(FindInNeighList(z, parent) >= 0);
                Unmerge(mi);
                EdgeNo eu = FindInNeighList(u, parent), ev = FindInNeighList(v, parent);
                Assert((eu >= 0 && ev < 0) || (eu < 0 && ev >= 0));
                if (eu < 0) { swap(eu, ev); swap(u, v); }
            }
            else
            {
                // Since z had no parent, this means that all edges incident on it pointed into z,
                // and they must now point into u or v.  It doesn't matter how we direct the edge
                // between u and v; let's say that u will become the parent of v.
                // - But, wait: since we do merges of degree-2 nodes before dealing with the degree-1 root,
                // it is possible that our current merged node is the result of merging the degree-1 root
                // with its (degree-2) neighbour; in this case we have to check if one of u and v has degree 1,
                // and we should make *it* the parent.
                Unmerge(mi);
                int uDeg = verts[u].deg, vDeg = verts[v].deg;
                Assert(uDeg > 1 || vDeg > 1);
                if (vDeg == 1) swap(u, v);
            }
            auto &U = verts[u], &V = verts[v]; U.parent = parent; V.parent = u;
            Assert(U.deg == 2 || V.deg == 2);
            for (EdgeNo e = U.firstEdge; e >= 0; ) { 
                Edge &E = edges[e]; VertNo w = E.Other(u); auto &W = verts[w];
                if (w == parent) { }
                else if (w == v) { Assert(W.parent == u); }
                else if (W.parent == z) W.parent = u; 
                e = E.Next(u); }
            for (EdgeNo e = V.firstEdge; e >= 0; ) { 
                Edge &E = edges[e]; VertNo w = E.Other(v); Assert(w != parent);
                auto &W = verts[w]; 
                if (w == u) { Assert(W.parent == parent); }
                else if (W.parent == z) W.parent = v;
                e = E.Next(v); }
        }
        // Now that we have unmerged everything, let's see if the leaves really
        // form a cycle and if the rest is really a tree.
        {
            // Find out which nodes are leaves; also find the root.
            vector<bool> isLeaf(verts.size(), true);
            VertNo root = -1;
            for (VertNo u = 0; u < (VertNo) verts.size(); ++u) { auto &U = verts[u]; if (! U.exists) continue;
                if (U.parent < 0) { Assert(root < 0); root = u; continue; }
                isLeaf[U.parent] = false; }
            Assert(root >= 0); Assert(! isLeaf[root]);
            // It's enough to check if every edge connects either two leaves or a node with its parent,
            // and that each leaf has degree 3.
            for (VertNo u = 0; u < (VertNo) verts.size(); ++u) { auto &U = verts[u]; if (! U.exists) continue;
                if (isLeaf[u]) { XAssert(U.deg == 3); }
                for (EdgeNo e = U.firstEdge; e >= 0; ) { auto &E = edges[e]; Assert(E.exists);
                    VertNo v = E.Other(u); Assert(verts[v].exists);
                    XAssert(v == U.parent || u == verts[v].parent || (isLeaf[u] && isLeaf[v])); 
                    e = E.Next(u); } }
            return true;
        }
    }

public:

    void RecoverTree()
    {
        // Any nodes of degree 2 must be internal nodes and can for the time being
        // be merged with one of their two neighbours.  The edges incident on such nodes
        // are guaranteed to not belong to the leaf cycle, and we'll preserve this
        // property later when merging; it can help us when choosing the root at the end.
        vector<MergeInfo> miTwos;
        for (VertNo u = 0; u < (VertNo) verts.size(); ++u) {
            auto &U = verts[u]; if (! U.exists) continue;
            Assert(U.deg >= 2); 
            if (U.deg != 2) continue;
            VertNo uNeigh = edges[U.firstEdge].Other(u);
            miTwos.push_back({}); auto &mi = miTwos.back();
            mi.inVerts.emplace_back(u); mi.inVerts.emplace_back(uNeigh); 
            Merge(mi); }
        // Now keep merging until we get rid of all the triangles.
        //const int vertsSizeBeforeTriangleMerging = verts.size();
        vector<MergeInfo> miThrees;
        for (VertNo u = 0; u < (VertNo) verts.size(); ++u)
        {
            auto &U = verts[u]; if (! U.exists) continue;
            Assert(U.deg >= 3);
            if (U.deg != 3) continue;
            // Prepare a list of u's three neighbours.
            VertNo uNeigh[3]; { EdgeNo e = U.firstEdge;
            for (int i = 0; i < 3; ++i) { 
                auto &E = edges[e]; Assert(E.exists); Assert(E.u == u || E.v == u);
                uNeigh[i] = E.Other(u); e = E.Next(u); }
            Assert(e < 0); }
            // For every two neighbours, say v and w, see if (u, v, w) form a triangle.
            for (int i = 0; i < 3; ++i) 
            {
                VertNo v = uNeigh[i], w = uNeigh[(i + 1) % 3];
                EdgeNo e = GetEdgeIfExists(v, w); if (e < 0) continue;
                // We have a triangle (u, v, w).  In any triangle, two of the three nodes
                // are leaves, the third is an internal node (except if there are only three
                // leaves left in the whole graph, in which case the leaf cycle is also a
                // triangle and consists of three leaves).  This means that at least two
                // of the three nodes u, v, w must have degree 3.  Since u already has it,
                // one of v and w must have it too.
                Assert(verts[v].exists); Assert(verts[w].exists);
                int vDeg = verts[v].deg, wDeg = verts[w].deg;
                Assert(vDeg == 3 || wDeg == 3);
                if (wDeg > 3) { swap(v, w); swap(vDeg, wDeg); }
                Assert(wDeg == 3);
                //                                                                    |                 |
                // It is possible that this triangle is one of a clump of             v                 v
                // two or more adjacent ones, as shown by the figure on             / | \      or     / | \         .
                // the right.  In this case one of v must have degree > 3,         u--w--z           w--u--y
                // but w must have degree exactly 3, and must have a third neighbour z which forms a triangle with v and w;
                // or u must have a third neighbour y which forms a triangle with v and u.  Note that we won't necessarily
                // note this second scenario when dealing with y because possibly u is the result of a recent merge
                // and did not yet exist when we were dealing with y.
                auto /*&V = verts[v],*/ &W = verts[w];
                VertNo z = -1, y = uNeigh[(i + 2) % 3]; int wDeg_ = 0; Assert(EdgeExists(u, y));
                for (EdgeNo e = W.firstEdge; e >= 0; ) {
                    auto &E = edges[e]; Assert(E.exists); ++wDeg_;
                    VertNo zz = E.Other(w); if (zz != u && zz != v) { Assert(z < 0); z = zz; }
                    e = E.Next(w); }
                Assert(wDeg_ == 3); Assert(z >= 0);
                if (vDeg > 3)
                {
                    bool clump = false; VertNo third = -1;
                    // Check if the triangle (v, w, z) exists.
                    // In that case we have a clump of triangles.  We'll just merge u and w for now.
                    // However, if an edge from u to z also exists, this means that u, w, z
                    // are the last three leaves, children of the root v - or in fact, since
                    // the situation is very symetrical (the graph is a tetrahedron at that point),
                    // it's hard to say what exactly the root is.  At any rate this would mean 
                    // that the merging is now at an end.  -  But wait, this shouldn't be possible
                    // since we asserted that v has degree > 3, so it can't be the root and there
                    // must be some other leaves elsewhere.
                    if (EdgeExists(z, v)) { clump = true; third = z; Assert( ! EdgeExists(u, z)); }
                    // Similarly check if the triangle (u, y, v) exists.
                    else if (EdgeExists(y, v)) { clump = true; third = y; Assert(! EdgeExists(w, y)); }
                    if (clump)
                    {
                        // Otherwise we'll merge u and w.  The merged node has two reasons to 
                        // be connected to v (edges u-v and v-w; or u-v and v-y), but our merging function creates just one edge.
                        // - By the way, at this point we also know that v must be the parent of u and w (or u and y),
                        // so the edges u-v and w-v (or y-v) aren't part of the leaf cycle.
                        { EdgeNo e = GetEdge(u, v); Assert(e >= 0); Assert(edges[e].exists); }
                        { EdgeNo e = GetEdge(third, v); Assert(e >= 0); Assert(edges[e].exists); }
                        miThrees.push_back({}); auto &mi = miThrees.back();
                        mi.inVerts.emplace_back(u); mi.inVerts.emplace_back(w); 
                        Merge(mi); break;
                    }
                }
                else 
                {
                    // u, v and w form a triangle and all have degree 3; hence we can merge the whole triangle.
                    // The only problem with this plan is if the whole graph has been reduced to a 
                    // tetrahedron by now, in which case z (w's third neighbour) is also u's and v's third neighbour.
                    // In that case we'll stop the merging.
                    if (EdgeExists(z, u) && EdgeExists(z, v)) { u = verts.size() + 1; break; }
                    miThrees.push_back({}); auto &mi = miThrees.back();
                    mi.inVerts.emplace_back(u); mi.inVerts.emplace_back(v); mi.inVerts.emplace_back(w);
                    Merge(mi); break;
                }
            } // for i
        } // for u
        // At this point what's left of the graph must be a tetrahedron, representing one
        // internal node with three leaves as children.
        {
            Assert(nVerts == 4);
            VertNo tVerts[4]; int nVerts_ = 0;
            for (VertNo u = 0; u < (VertNo) verts.size(); ++u) { 
                auto &U = verts[u]; if (! U.exists) continue;
                tVerts[nVerts_++] = u; 
                Assert(U.deg == 3); }
            Assert(nVerts_ == 4);
            for (int i = 1; i < 4; ++i) for (int j = 0; j < i; ++j) Assert(EdgeExists(tVerts[j], tVerts[i]));
            for (int i = 0; i < 4; ++i)
            {
                auto backup = *this;
                if (RecoverTree_TryRoot(tVerts[i], tVerts, miThrees, miTwos)) return;
                *this = backup;
            }
            Assert(false);
        } // unmerging
    } // RecoverTree
};

enum { MaxN = 100'000, MaxK = 100'000, MaxKSmall = 10, M = 1'000'000'007 };

typedef enum { cc000, cc001, cc010, cc100, cc012  } ColorCombo;
inline ColorCombo Normalize(int b1, int b2, int b3)
{
    if (b1 == b2) return (b3 == b1) ? cc000 : cc001;
    else if (b2 == b3) return cc100;
    else return (b1 == b3) ? cc010 : cc012;
}

struct Solution2
{
    int f[5] = {}, g[5] = {};
    int &F(int b1, int b2, int b3) { return f[Normalize(b1, b2, b3)]; }
    int &G(int b1, int b2, int b3) { return g[Normalize(b1, b2, b3)]; }
    int F(int b1, int b2, int b3) const { return f[Normalize(b1, b2, b3)]; }
    int G(int b1, int b2, int b3) const { return g[Normalize(b1, b2, b3)]; }
};

// Assuming that a and b are in the range 0..M-1, this function returns (a - b) mod M.
inline int SubM(int a, int b) { if (b > 0) { a += (M - b); if (a >= M) a -= M; } return a; }
inline llint MulM(llint a, llint b) { return ((a % M) * (b % M)) % M; }

struct Vert2
{
    VertNo parent = -1, firstChild = -1, nextSib = -1;
    VertNo firstLeaf = -1, lastLeaf = -1;
    Solution2 sol2;
};

struct Graph2
{
    int n = 0, k = -1; VertNo root = -1;
    vector<Vert2> verts;
    vector<VertNo> prevLeaf, nextLeaf;

private:

    void GetPreorder(vector<VertNo> &preorder)
    {
        preorder.clear(); preorder.emplace_back(root);
        for (int i = 0; i < (int) preorder.size(); ++i) {
            VertNo u = preorder[i]; auto &U = verts[u];
            for (VertNo v = U.firstChild; v >= 0; v = verts[v].nextSib)
                preorder.emplace_back(v); }
        Assert(int(preorder.size()) == n);
    }

public:

    void Init(const Graph &G)
    {
        // Prepare the vector of vertices with provisional (unordered) child lists.
        // Also prepare a vector indicating which vertices are leaves.
        n = G.nVerts; k = -1; verts.clear(); verts.resize(n);
        vector<bool> isLeaf(n, true); root = -1;
        for (VertNo u = 0; u < (VertNo) G.verts.size(); ++u) 
        {
            auto &GU = G.verts[u]; Assert(u < n ? GU.exists : ! GU.exists);
            if (u >= n) continue;
            auto &U = verts[u]; U.parent = GU.parent;
            if (GU.parent < 0) { Assert(root < 0); root = u; continue; }
            isLeaf[U.parent] = false;
            auto &P = verts[U.parent]; U.nextSib = P.firstChild; P.firstChild = u;
        }
        Assert(root >= 0);
        int nLeaves = 0; for (VertNo u = 0; u < n; ++u) if (isLeaf[u]) ++nLeaves;
        // Prepare a preorder traversal.
        vector<VertNo> preorder; GetPreorder(preorder);
        // For each leaf, determine the previous and next leaf on the leaf cycle.
        prevLeaf.clear(); prevLeaf.resize(n, -1);
        nextLeaf.clear(); nextLeaf.resize(n, -1);
        VertNo firstLeaf = 0; while (firstLeaf < n && ! isLeaf[firstLeaf]) ++firstLeaf;
        Assert(firstLeaf < n); VertNo lastLeaf = -1;
        int nLeaves2 = 0;
        for (VertNo u = firstLeaf, prev = -1; ; )
        {
            auto &GU = G.verts[u]; Assert(GU.deg == 3);
            auto &U = verts[u]; VertNo parent = U.parent; Assert(parent >= 0);
            ++nLeaves2;
            bool seenPrev = false, seenParent = false; VertNo next = -1;
            for (EdgeNo e = GU.firstEdge; e >= 0; ) { auto &E = G.edges[e];
                VertNo v = E.Other(u); 
                if (v == parent) { Assert(! seenParent); seenParent = true; }
                else if (v == prev) { Assert(! seenPrev); seenPrev = true; }
                else if (next < 0) next = v;
                else { Assert(u == firstLeaf); lastLeaf = v; }
                e = E.Next(u); }
            Assert(seenParent); if (prev >= 0) Assert(seenPrev); 
            Assert(next >= 0); 
            if (next == firstLeaf) { Assert(u == lastLeaf); }
            if (u == firstLeaf) { Assert(lastLeaf >= 0); prevLeaf[u] = lastLeaf; }
            else { Assert(prev >= 0); prevLeaf[u] = prev; }
            nextLeaf[u] = next; prev = u; u = next;
            if (u == firstLeaf) break;
        }
        Assert(nLeaves == nLeaves2);
        // Prepare the ordered lists of children.
        vector<VertNo> temp(n, -1); vector<int> nLeavesInSubtree(n, 0);
        for (int i = int(preorder.size()) - 1; i >= 0; --i)
        {
            VertNo u = preorder[i]; auto &U = verts[u];
            vector<VertNo> children; for (VertNo v = U.firstChild; v >= 0; v = verts[v].nextSib) children.emplace_back(v);
            if (children.empty()) { U.firstLeaf = u; U.lastLeaf = u; nLeavesInSubtree[u] = 1; continue; }
            // For each child 'v', let temp[v.firstLeaf] and temp[v.lastLeaf] point to 'v'.
            for (VertNo v : children) { auto &V = verts[v]; Assert(V.parent == u);
                V.nextSib = -1; nLeavesInSubtree[u] += nLeavesInSubtree[v];
                temp[V.firstLeaf] = -1; temp[prevLeaf[V.firstLeaf]] = -1;
                temp[V.lastLeaf] = -1; temp[nextLeaf[V.lastLeaf]] = -1; }
            for (VertNo v : children) { auto &V = verts[v];
                Assert(temp[V.firstLeaf] < 0); Assert(temp[V.lastLeaf] < 0);
                temp[V.firstLeaf] = v; temp[V.lastLeaf] = v; }
            // For each child v, we can use temp[prevLeaf[v.firstLeaf]] to determine the previous child
            // and temp[nextLeaf[v.lastLeaf]] to determine the next child.
            U.firstChild = -1; U.firstLeaf = -1; U.lastLeaf = -1;
            for (VertNo v : children) { auto &V = verts[v];
                VertNo prev = temp[prevLeaf[V.firstLeaf]];
                if (prev < 0) { Assert(U.firstChild < 0); U.firstChild = v; U.firstLeaf = V.firstLeaf; }
                V.nextSib = temp[nextLeaf[V.lastLeaf]];
                if (V.nextSib < 0) { Assert(U.lastLeaf < 0); U.lastLeaf = V.lastLeaf; } }
            if (U.firstChild >= 0) { Assert(U.firstLeaf >= 0); Assert(U.lastLeaf >= 0); Assert(nLeavesInSubtree[u] < nLeaves); }
            else {
                // If our children cover all the leaves, the previous for loop has actually
                // linked them into a cycle; we'll now interrupt the cycle at an arbitrary point.
                Assert(U.firstLeaf < 0); Assert(U.lastLeaf < 0); Assert(nLeavesInSubtree[u] == nLeaves);
                VertNo v = children[0]; auto &V = verts[v]; U.firstLeaf = V.firstLeaf; U.firstChild = v;
                VertNo w = temp[prevLeaf[V.firstLeaf]]; Assert(w >= 0);
                auto &W = verts[w]; Assert(W.nextSib == v); W.nextSib = -1;
                if (children.size() == 1) Assert(w == v); 
                U.lastLeaf = W.lastLeaf;
            }
        }
    }

    void CalculateVertex2(VertNo u) // O(1) time and space
    {
        auto &U = verts[u]; 
        auto &sU = U.sol2; 
        thread_local const Solution2 dummySol;
        VertNo v = U.firstChild; const Solution2 &sV = (v >= 0) ? verts[v].sol2 : dummySol;
        VertNo w = U.nextSib; const Solution2 &sW = (w >= 0) ? verts[w].sol2 : dummySol;
        // Let T(u) be the subgraph induced by u and all of its descendants;
        // and let T'(u) be the subgraph induced by u, all of its right siblings,
        // and all the descendants of u and these siblings.
        //   f_u(b1, b2, b3) = the number of colorings of T(u) where u has color b1,
        //                     the leftmost leaf of T(u) has color b2 and the rightmost leaf has color b3.
        //   g_u(b1, b2, b3) = the number of colorings of T'(u) where u's parent is presumed to have color b1
        //                     (this parent isn't a part of T'(u)), meaning that u and its siblings can't have color b1;
        //                     and where the leftmost leaf of T'(u) has color b2 and the rightmost leaf has color b3.
        // Edge case for f_u: if u is a leaf, then b1, b2, b3 all refer to node u, so there is 1 coloring
        // if b1 == b2 == b3, otherwise none.  
        // Otherwise, if u has the first child v, we have f_u(b1, b2, b3) = g_v(b1, b2, b3).
        if (v < 0) {
            //for (int b1 = 0; b1 < k; ++b1) for (int b2 = 0; b2 < k; ++b2) for (int b3 = 0; b3 < k; ++b3) fu[b1][b2][b3] = (b1 == b2 && b2 == b3) ? 1 : 0;
            for (int q = 0; q < 5; ++q) sU.f[q] = 0;
            sU.f[cc000] = 1;
        }
        else {
            //auto &gv = Vsol->g;
            //for (int b1 = 0; b1 < k; ++b1) for (int b2 = 0; b2 < k; ++b2) for (int b3 = 0; b3 < k; ++b3) fu[b1][b2][b3] = gv[b1][b2][b3]; }
            for (int q = 0; q < 5; ++q) sU.f[q] = sV.g[q]; }
        //
        // Now let's consider g_u.  We'll denote u's right sibling by w; let x be the rightmost leaf of T(u)
        // and let y be the leftmost leaf of T'(w).  Let bu, bx, and by be the color of u, x, and y, respectively.
        // Then bx and by must be different, because x and y are neighbours on the leaf cycle; and bu must be different from b1
        // by the definition of g_u.
        // We then have f_u(bu, b2, bx) ways of colouring T(u) and g_w(b1, by, b3) ways of colouring T'(w).
        // - An edge case: if u is the rightmost sibling, there is no w.
        //   Then g_u(b1, b2, b3) = sum_{bu} [bu != b1] f_u(bu, b2, b3).
        //thread_local int sum1[MaxK][MaxK], sum2[MaxK][MaxK];
        if (w < 0) {
            // sum1[0][0] = sum_{bu} f_u(bu, 0, 0) = f_u(0, 0, 0) + (k - 1) f_u(1, 0, 0).
            int sum1_00 = int((llint(sU.f[cc000]) + (k - 1) * llint(sU.f[cc100])) % M);
            // sum1[0][1] = sum_{bu} f_u(bu, 0, 1) = f_u(0, 0, 1) + f_u(1, 0, 1) + (k - 2) f_u(2, 0, 1).
            int sum1_01 = int((llint(sU.f[cc001]) + llint(sU.f[cc010]) + (k - 2) * llint(sU.f[cc012])) % M);
            // Now sum1[b2][b3] = sum_{bu} f_u(bu, b2, b3).
            // g_u(b1, b2, b3) = sum1(b2, b3) - f_u(b1, b2, b3)
            sU.g[cc000] = SubM(sum1_00, sU.f[cc000]);
            if (k > 1) {
                sU.g[cc001] = SubM(sum1_01, sU.f[cc001]);
                sU.g[cc010] = SubM(sum1_01, sU.f[cc010]);
                sU.g[cc100] = SubM(sum1_00, sU.f[cc100]);
                if (k > 2) sU.g[cc012] = SubM(sum1_01, sU.f[cc012]); }
        }
        else if (v < 0) {
            // u has a right sibling, but u itself is a leaf.  This means that b2, bx and bu all refer to the same node.
            //     g_u(b1, b2, b3) = sum_{bu, bx, by} [bu != b1, bx != by] f_u(bu, b2, bx) g_w(b1, by, b3)         // the general case
            //                     = sum_{by} [b2 != b1, b2 != by] f_u(b2, b2, b2) g_w(b1, by, b3)             // require bu = bx = b2 because they refer to the same node, namely u
            //                     = [b2 != b1] sum_{by} [b2 != by] g_w(b1, by, b3)  
            // sum1[0][0] = sum_{by} g_w(0, by, 0) = g_w(0, 0, 0) + (k - 1) g_w(0, 1, 0).
            int sum1_00 = int((llint(sW.g[cc000]) + (k - 1) * llint(sW.g[cc010])) % M);
            // sum1[0][1] = sum_{by} g_w(0, by, 1) = g_w(0, 0, 1) + g_w(0, 1, 1) + (k - 2) g_w(0, 2, 1).
            int sum1_01 = int((llint(sW.g[cc001]) + llint(sW.g[cc100]) + (k - 2) * llint(sW.g[cc012])) % M);
            // Now sum1[b1][b2] = sum_{by} g_w(b1, by, b3).
            // g_u(b1, b2, b3) = [b2 != b1] (sum1(b1, b3) - g_w(b1, b2, b3)).
            sU.g[cc000] = 0; // b1 == b2
            sU.g[cc001] = 0; // b1 == b2
            sU.g[cc010] = SubM(sum1_00, sW.g[cc010]);
            sU.g[cc100] = SubM(sum1_01, sW.g[cc100]);
            sU.g[cc012] = SubM(sum1_01, sW.g[cc012]);
        }
        else
        {
            // Let's consider the general case.  We have
            //     g_u(b1, b2, b3) = sum_{bu, bx, by} [bu != b1, bx != by] f_u(bu, b2, bx) g_w(b1, by, b3)
            //                     = sum_{bx} (sum_bu [bu != b1] f_u(bu, b2, bx) ) ( sum_by [bx != by]  g_w(b1, by, b3) )
            // So it will be useful to define sum1[b2][bx] = sum_bu f_u(bu, b2, bx)
            //                            and sum2[b1][b3] = sum_by g_w(b1, by, b3).
            // sum1(0, 0) = sum_bu f_u(bu, 0, 0) = f_u(0, 0, 0) + (k - 1) f_u(1, 0, 0).
            int sum1_00 = int((llint(sU.f[cc000]) + (k - 1) * llint(sU.f[cc100])) % M);
            // sum1(0, 1) = sum_bu f_u(bu, 0, 1) = f_u(0, 0, 1) + f_u(1, 0, 1) + (k - 2) f_u(2, 0, 1).
            int sum1_01 = int((llint(sU.f[cc001]) + llint(sU.f[cc010]) + (k - 2) * llint(sU.f[cc012])) % M);
            // sum2(0, 0) = sum_by g_w(0, by, 0) = g_w(0, 0, 0) + (k - 1) g_w(0, 1, 0).
            int sum2_00 = int((llint(sW.g[cc000]) + (k - 1) * llint(sW.g[cc010])) % M);
            // sum2(0, 1) = sum_by g_w(0, by, 1) = g_w(0, 0, 1) + g_w(0, 1, 1) + (k - 2) g_w(0, 2, 1).
            int sum2_01 = int((llint(sW.g[cc001]) + llint(sW.g[cc100]) + (k - 2) * llint(sW.g[cc012])) % M);
            for (int cc = 0; cc < 5; ++cc) {
                int b1 = (cc == cc100) ? 1 : 0, b2 = (cc == cc010 || cc == cc012) ? 1 : 0, b3 = (cc == cc012) ? 2 : (cc == cc001) ? 1 : 0;
                //  g_u(b1, b2, b3) = sum_bx [sum1(b2, bx) - f_u(b1, b2, bx)] * [sum2(b1, b3) - g_w(b1, bx, b3)].
                llint result = 0;
                for (int bx = 0; bx < k && bx <= 3; ++bx) {
                    llint term1 = SubM((bx == b2 ? sum1_00 : sum1_01), sU.F(b1, b2, bx));
                    llint term2 = SubM((b1 == b3 ? sum2_00 : sum2_01), sW.G(b1, bx, b3));
                    // When bx >= 3, it is certainly different from b1, b2 and b3; therefore the
                    // values of term1 and term2 will be the same for bx = 3, 4, ..., k-1,
                    // so we can simply take the product for bx = 3 and multiply it by k-3.
                    llint product = (term1 * term2) % M; if (bx == 3 && k > 4) product = (product * (k - 3)) % M;
                    result += product; }
                sU.g[cc] = int(result % M);
            }
        }
        // If u and its right siblings cover all the leaves, then the leaves referred to by b2 and b3 are connected.
        // Thus if b2 and b3 are equal, no valid coloring is possible.
        VertNo lastLeaf = (U.parent >= 0) ? verts[U.parent].lastLeaf : U.lastLeaf;
        if (prevLeaf[U.firstLeaf] == lastLeaf) 
            //for (int b1 = 0; b1 < k; ++b1) for (int b2 = 0; b2 < k; ++b2) gu[b1][b2][b2] = 0;
            sU.g[cc000] = 0, sU.g[cc100] = 0;
    }

public:

    int Calculate(int k_)
    {
        k = k_; Assert(1 <= k); Assert(k <= MaxK);
        vector<VertNo> preorder; GetPreorder(preorder);
        for (int i = int(preorder.size()) - 1; i >= 0; --i) {
            VertNo u = preorder[i]; auto &U = verts[u];
            CalculateVertex2(u); }
        Assert(root >= 0); auto &R = verts[root]; 
        auto &sR = R.sol2;
        // We want to return sum_b1 sum_b2 sum_b3 f_root(b1, b2, b3).
        // 000: can be chosen in k ways (choose b1 from 0 to k-1, then set b2 = b1 and b3 = b1);
        // 001, 010, 100: can be chosen in k * (k - 1) ways
        // 012: can be chosen in k * (k - 1) * (k - 2) ways
        int total2 = int(MulM(MulM(MulM(sR.f[cc012], k - 2) + llint(sR.f[cc001]) + llint(sR.f[cc010]) + llint(sR.f[cc100]), k - 1) + sR.f[cc000], k));
        return total2;
    }
};

int main()
{
    // Read the input graph.
    int n, m, k; int ok_ = scanf("%d %d %d", &n, &m, &k); Assert(ok_ == 3);
    Assert(3 <= n); Assert(n <= MaxN); Assert(3 <= m); 
    Assert(m <= 2 * (n - 1)); // the tree had n-1 edges and at most n-1 leaves, therefore the leaf cycle had at most n-1 edges, for a total of at most 2(n-1) edges
    Assert(1 <= k); Assert(k <= MaxK);
    Assert(m >= n + 2); // n-1 edges for the tree + as least 3 more for the leaf cycle
    Graph G; for (VertNo u = 0; u < n; ++u) { VertNo uu = G.AddVert(); Assert(uu == u); }
    for (EdgeNo e = 0; e < m; ++e) {
        int u, v; ok_ = scanf("%d %d", &u, &v); Assert(ok_ == 2);
        Assert(1 <= u); Assert(u <= n); Assert(1 <= v); Assert(v <= n); Assert(u != v);
        --u; --v; EdgeNo ee = G.AddEdge(u, v); Assert(ee == e); }
    // Divide it into the tree and the leaf cycle.
    G.RecoverTree();
    // Convert the tree into an ordered tree.
    Graph2 G2; G2.Init(G);
    // Calculate the number of colourings and output it.
    printf("%d\n", G2.Calculate(k));
    return 0;
}