/**
 * @file        graph.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#include <dfy/graph.hpp>

#include <set>
#include <queue>


dfy::Graph::Graph(const std::vector<std::pair<int, int>>& Edges,
                  const std::vector<double>& Weights)
{
    // Get sorted edges and number of nodes
    std::map<std::pair<int, int>, double> EdgeSet;
    int nNodes = 0;
    for (int i = 0; i < Edges.size(); ++i)
    {
        std::pair<int, int> e = Edges[i];
        // We don't care about self-edges
        if (e.first == e.second)
        {
            nNodes = std::max(nNodes, e.first + 1);
            continue;
        }
        double w = 1.0;
        if (!Weights.empty())
            w = Weights[i];
        EdgeSet.emplace(e, w);
        std::swap(e.first, e.second);
        EdgeSet.emplace(e, w);

        nNodes = std::max(nNodes, std::max(e.first, e.second) + 1);
    }

    // Initialize indices and allocate edge array
    m_Idxs.resize(nNodes + 1);
    m_Adjs.reserve(EdgeSet.size());

    // Populate adjacency list
    int CurNode = 0;
    for (auto ew : EdgeSet)
    {
        std::pair<int, int> e = ew.first;
        double w = ew.second;
        // If node has changed, we update offsets until
        // we reach the current node
        while (CurNode != e.first)
            m_Idxs[++CurNode] = m_Adjs.size();
        m_Adjs.emplace_back(e.second, w);
    }
    for (CurNode++; CurNode <= nNodes; ++CurNode)
        m_Idxs[CurNode] = m_Adjs.size();
}

dfy::Graph::Graph(const dfy::Graph& G)
{
    m_Idxs = G.m_Idxs;
    m_Adjs = G.m_Adjs;
}

dfy::Graph::Graph(dfy::Graph&& G)
{
    m_Idxs = std::move(G.m_Idxs);
    m_Adjs = std::move(G.m_Adjs);
}

dfy::Graph& dfy::Graph::operator=(const dfy::Graph& G)
{
    m_Idxs = G.m_Idxs;
    m_Adjs = G.m_Adjs;
    return *this;
}

dfy::Graph& dfy::Graph::operator=(dfy::Graph&& G)
{
    m_Idxs = std::move(G.m_Idxs);
    m_Adjs = std::move(G.m_Adjs);
    return *this;
}

dfy::Graph::~Graph() { }

int dfy::Graph::NumNodes() const { return m_Idxs.size() - 1; }
int dfy::Graph::NumAdjacents(int i) const { return m_Idxs[i + 1] - m_Idxs[i]; }
int dfy::Graph::NumEdges() const { return m_Adjs.size(); }

const dfy::WEdge &dfy::Graph::GetAdjacent(int node_i, int adj_i) const
{
    return m_Adjs[m_Idxs[node_i] + adj_i];
}

dfy::Path dfy::Graph::DijkstraPath(int src, int dst) const
{
    Eigen::VectorXd Dists;
    Dists.setConstant(NumNodes(), std::numeric_limits<double>::infinity());

    Eigen::VectorXi Parent;
    Parent.setConstant(NumNodes(), -1);

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    Dists[src] = 0.0;
    Q.emplace(0.0, src);
    while (!Q.empty())
    {
        int i;
        double wi;
        std::tie(wi, i) = Q.top();
        Q.pop();

        if (i == dst)
            break;
        
        int Degree = NumAdjacents(i);
        for (int jj = 0; jj < Degree; ++jj)
        {
            int j;
            double wj;
            std::tie(j, wj) = GetAdjacent(i, jj);
            if (Dists[j] <= wi + wj)
                continue;
            Dists[j] = wi + wj;
            Parent[j] = i;
            Q.emplace(Dists[j], j);
        }
    }

    std::vector<dfy::WEdge> Path;
    double Length = Dists[dst];
    while (Parent[dst] != -1)
    {
        Path.emplace_back(dst, Dists[dst] - Dists[Parent[dst]]);
        dst = Parent[dst];
    }
    Path.emplace_back(dst, Dists[dst]);
    std::reverse(Path.begin(), Path.end());
    return { Length, Path };
}

Eigen::VectorXd dfy::Graph::DijkstraDistance(int src) const
{
    Eigen::VectorXd D;
    DijkstraDistance(src, D);
    return D;
}

void dfy::Graph::DijkstraDistance(int src, Eigen::VectorXd& Dists) const
{
    Dists.setConstant(NumNodes(), std::numeric_limits<double>::infinity());

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    Dists[src] = 0.0;
    Q.emplace(0.0, src);
    while (!Q.empty())
    {
        int i = Q.top().second;
        double wi = Dists[i];
        // std::tie(wi, i) = Q.top();

        Q.pop();
        
        int Degree = NumAdjacents(i);
        for (int jj = 0; jj < Degree; ++jj)
        {
            int j;
            double wj;
            std::tie(j, wj) = GetAdjacent(i, jj);
            if (Dists[j] <= wi + wj)
                continue;
            Dists[j] = wi + wj;
            Q.emplace(Dists[j], j);
        }
    }
}

Eigen::VectorXd dfy::Graph::DijkstraDistance(int src, const Eigen::VectorXi &Tag, int Filter) const
{
    Eigen::VectorXd Dists;
    Dists.setConstant(NumNodes(), std::numeric_limits<double>::infinity());

    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    Dists[src] = 0.0;
    Q.emplace(0.0, src);
    while (!Q.empty())
    {
        int i;
        double wi;
        std::tie(wi, i) = Q.top();
        Q.pop();
        
        int Degree = NumAdjacents(i);
        for (int jj = 0; jj < Degree; ++jj)
        {
            int j;
            double wj;
            std::tie(j, wj) = GetAdjacent(i, jj);
            if (Tag[j] != Filter)
                continue;
            if (Dists[j] <= wi + wj)
                continue;
            Dists[j] = wi + wj;
            Q.emplace(Dists[j], j);
        }
    }

    return Dists;
}

std::vector<int> dfy::Graph::ConnectedComponents() const
{
    std::vector<int> CC;
    CC.resize(NumNodes(), 0);
    int CurCC = 0;

    Eigen::VectorXi Visited;
    Visited.resize(NumNodes());
    Visited.setZero();

    while (Visited.sum() < NumNodes())
    {
        int Root = 0;
        for (int i = 0; i < NumNodes(); ++i)
        {
            if (Visited[i] == 0)
            {
                Root = i;
                break;
            }
        }

        std::queue<int> Q;
        Q.emplace(Root);
        while (!Q.empty())
        {
            int N = Q.front();
            Q.pop();
            
            if (Visited[N] != 0)
                continue;

            CC[N] = CurCC;
            Visited[N] = 1;
            int DegN = NumAdjacents(N);
            for (int jj = 0; jj < DegN; ++jj)
            {
                int j = GetAdjacent(N, jj).first;
                Q.emplace(j);
            }
        }

        CurCC += 1;
    }

    return CC;
}

std::vector<std::pair<int, int>> dfy::Graph::MinSpanTree() const
{
    std::vector<std::pair<int, int>> MST;
    MST.reserve(NumNodes() - 1);

    std::vector<bool> Visited;
    Visited.resize(NumNodes());

    std::priority_queue<std::pair<double, std::pair<int, int>>> Q;
    std::pair<int, int> Start{ 0, -1 };
    Q.emplace(0.0, Start);
    while (!Q.empty())
    {
        double w;
        std::pair<int, int> v;
        std::tie(w, v) = Q.top();
        Q.pop();
        if (Visited[v.first])
            continue;
        Visited[v.first] = true;
        if (v.second != -1)
            MST.push_back(v);

        int nadj = NumAdjacents(v.first);
        for (int ii = 0; ii < nadj; ++ii)
        {
            int i;
            double wi;
            std::tie(i, wi) = GetAdjacent(v.first, ii);
            std::pair<int, int> vi{ i, v.first };
            Q.emplace(wi, vi);
        }
    }

    return MST;
}

dfy::Graph dfy::Graph::SubGraph(const std::vector<int> &Indices) const
{
    std::map<int, int> IdxInv;
    for (int i = 0; i < Indices.size(); ++i)
        IdxInv.emplace(Indices[i], i);
    
    std::vector<std::pair<int, int>> Edges;
    std::vector<double> Weights;
    Edges.reserve(NumEdges());
    Weights.reserve(NumEdges());

    for (auto it : IdxInv)
    {
        int iOld = it.first;
        int iNew = it.second;
        int nadj = NumAdjacents(iOld);
        // To ensure we don't fuck up disconnected components
        bool Isolated = true;
        for (int jj = 0; jj < nadj; ++jj)
        {
            int jOld;
            double w;
            std::tie(jOld, w) = GetAdjacent(iOld, jj);
            if (IdxInv.find(jOld) == IdxInv.end())
                continue;
            Isolated = false;
            int jNew = IdxInv[jOld];
            Edges.emplace_back(iNew, jNew);
            Weights.emplace_back(w);
        }
        if (Isolated)
        {
            Edges.emplace_back(iNew, iNew);
            Weights.emplace_back(std::numeric_limits<double>::infinity());
        }
    }

    return dfy::Graph(Edges, Weights);
}


dfy::Graph dfy::MeshToGraph(const dfy::Mesh &M, M2GDist Dist)
{
    std::vector<std::pair<int, int>> Edges;
    Edges.reserve(3 * M.NumTriangles());
    for (int i = 0; i < M.NumTriangles(); ++i)
    {
        std::pair<int, int> e;

        e.first = M.Triangles()(i, 0);
        e.second = M.Triangles()(i, 1);
        if (e.first > e.second)
            std::swap(e.first, e.second);
        Edges.emplace_back(e);
        
        e.first = M.Triangles()(i, 1);
        e.second = M.Triangles()(i, 2);
        if (e.first > e.second)
            std::swap(e.first, e.second);
        Edges.emplace_back(e);
        
        e.first = M.Triangles()(i, 2);
        e.second = M.Triangles()(i, 0);
        if (e.first > e.second)
            std::swap(e.first, e.second);
        Edges.emplace_back(e);
    }

    std::sort(Edges.begin(), Edges.end());
    auto EEnd = std::unique(Edges.begin(), Edges.end());
    Edges.erase(EEnd, Edges.end());

    std::vector<double> Weights;
    Weights.resize(Edges.size());
    for (int i = 0; i < Edges.size(); ++i)
        Weights[i] = Dist(M, Edges[i].first, Edges[i].second);
    
    return dfy::Graph(Edges, Weights);
}

dfy::Graph dfy::DualMeshToGraph(const dfy::ManifoldMesh &M, DM2GDist Dist)
{
    std::vector<std::pair<int, int>> Edges;
    Edges.reserve(3 * M.NumTriangles());
    for (int i = 0; i < M.NumTriangles(); ++i)
    {
        std::pair<int, int> e;

        e.first = i;
        e.second = M.TriTriAdj()(i, 0);
        if (e.second >= 0)
        {
            if (e.first > e.second)
                std::swap(e.first, e.second);
            Edges.emplace_back(e);
        }

        e.first = i;
        e.second = M.TriTriAdj()(i, 1);
        if (e.second >= 0)
        {
            if (e.first > e.second)
                std::swap(e.first, e.second);
            Edges.emplace_back(e);
        }

        e.first = i;
        e.second = M.TriTriAdj()(i, 2);
        if (e.second >= 0)
        {
            if (e.first > e.second)
                std::swap(e.first, e.second);
            Edges.emplace_back(e);
        }
    }

    std::sort(Edges.begin(), Edges.end());
    auto EEnd = std::unique(Edges.begin(), Edges.end());
    Edges.erase(EEnd, Edges.end());

    std::vector<double> Weights;
    Weights.resize(Edges.size());
    for (int i = 0; i < Edges.size(); ++i)
        Weights[i] = Dist(M, Edges[i].first, Edges[i].second);
    
    return dfy::Graph(Edges, Weights);
}

double dfy::EuclideanDistance(const dfy::Mesh &M, int i, int j)
{
    return (M.Vertices().row(i) - M.Vertices().row(j)).norm();
}

double dfy::DualEuclideanDistance(const dfy::ManifoldMesh &M, int i, int j)
{
    return (M.FaceBarycs().row(i) - M.FaceBarycs().row(j)).norm();
}

double dfy::AngularDistance(const dfy::Mesh &M, int i, int j)
{
    double d = std::acos(M.VertNormals().row(i).dot(M.VertNormals().row(j)));
    // return std::min(2.0, std::max(0.0, 1 - d));
    d = std::min(1.0, std::max(-1.0, d));
    return std::acos(d);
}

double dfy::DualAngularDistance(const dfy::ManifoldMesh &M, int i, int j)
{
    double d = M.FaceNormals().row(i).dot(M.FaceNormals().row(j));
    return std::min(2.0, std::max(0.0, 1 - d));
    // d = std::min(1.0, std::max(-1.0, d));
    // return std::acos(d);
}

double dfy::GeodesicDistance(const dfy::ManifoldMesh &M, int i, int j)
{
    int i1, i2, i3;
    for (i1 = 0; i1 < 3; ++i1)
    {
        if (M.TriTriAdj()(i, i1) == j)
        {
            i2 = (i1 + 1) % 3;
            i3 = (i1 + 2) % 3;
        }
    }
    int j1, j2, j3;
    for (j1 = 0; j1 < 3; ++j1)
    {
        if (M.TriTriAdj()(j, j1) == i)
        {
            j2 = (j1 + 1) % 3;
            j3 = (j1 + 2) % 3;
        }
    }
    //   i3-j2
    //  /  |  \
    // i1  |  j1
    //  \  |  /
    //   i2|j3
    double theta = M.FaceAngles()(i, i2) + M.FaceAngles()(j, j3);
    double a = M.EdgeLengths()(i, i3);
    double b = M.EdgeLengths()(j, j2);
    // bi = (i1 + i2 + i3) / 3;
    // bj = (j1 + j2 + j3) / 3 = (j1 + i3 + i2) / 3
    // bi - bj = (i1 - j1) / 3
    // || bi - bj || = || i1 - j1 || / 3
    double d = std::sqrt(a * a + b * b - 2 * a * b * std::cos(theta)) / 3.0;
    if (std::isnan(d))
        return dfy::DualEuclideanDistance(M, i, j);
    return d;
}

double dfy::ConstantDistance(const dfy::Mesh &, int i, int j)
{
    return 1.0;
}

double dfy::DualConstantDistance(const dfy::ManifoldMesh &, int i, int j)
{
    return 1.0;
}