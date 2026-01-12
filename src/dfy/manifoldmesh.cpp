/**
 * @file        manifoldmesh.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#include <dfy/manifoldmesh.hpp>
#include <dfy/utils.hpp>


dfy::ManifoldMesh::ManifoldMesh(const Eigen::MatrixXd& Verts,
                                const Eigen::MatrixXi& Tris)
    : dfy::Mesh(Verts, Tris)
{
    InitManifoldMesh();
}

dfy::ManifoldMesh::ManifoldMesh(const std::string& Filename)
    : dfy::Mesh(Filename)
{
    InitManifoldMesh();
}

dfy::ManifoldMesh::ManifoldMesh(const dfy::ManifoldMesh& M)
    : dfy::Mesh(M)
{
    m_E = M.m_E;
    m_EL = M.m_EL;
    m_E2T = M.m_E2T;
    m_T2E = M.m_T2E;
    m_T2T = M.m_T2T;
}

dfy::ManifoldMesh::ManifoldMesh(dfy::ManifoldMesh&& M)
    : dfy::Mesh(M)
{
    m_E = std::move(M.m_E);
    m_EL = std::move(M.m_EL);
    m_E2T = std::move(M.m_E2T);
    m_T2E = std::move(M.m_T2E);
    m_T2T = std::move(M.m_T2T);
}

dfy::ManifoldMesh& dfy::ManifoldMesh::operator=(const dfy::ManifoldMesh& M)
{
    dfy::Mesh::operator=(M);
    m_E = M.m_E;
    m_EL = M.m_EL;
    m_E2T = M.m_E2T;
    m_T2E = M.m_T2E;
    m_T2T = M.m_T2T;
    return *this;
}

dfy::ManifoldMesh& dfy::ManifoldMesh::operator=(dfy::ManifoldMesh&& M)
{
    dfy::Mesh::operator=(M);
    m_E = std::move(M.m_E);
    m_EL = std::move(M.m_EL);
    m_E2T = std::move(M.m_E2T);
    m_T2E = std::move(M.m_T2E);
    m_T2T = std::move(M.m_T2T);
    return *this;
}

dfy::ManifoldMesh::~ManifoldMesh() { }

int dfy::ManifoldMesh::NumEdges() const { return m_E.rows(); }
const Eigen::MatrixXi &dfy::ManifoldMesh::Edges() const { return m_E; }
const Eigen::VectorXd &dfy::ManifoldMesh::EdgeLengths() const { return m_EL; }
const Eigen::MatrixXi &dfy::ManifoldMesh::EdgeTriAdj() const { return m_E2T; }
const Eigen::MatrixXi &dfy::ManifoldMesh::TriEdgeAdj() const { return m_T2E; }
const Eigen::MatrixXi &dfy::ManifoldMesh::TriTriAdj() const { return m_T2T; }

Eigen::MatrixXd dfy::ManifoldMesh::EdgeCenters() const
{
    Eigen::MatrixXd EC = Vertices()(m_E(Eigen::all, 0), Eigen::all);
    EC += Vertices()(m_E(Eigen::all, 1), Eigen::all);
    EC = 0.5 * EC.array();
    return EC;
}



void dfy::ManifoldMesh::InitManifoldMesh()
{
    // Create vector of edges
    // Each edge keeps track of its generating triangle and 
    // its original index in the edges array
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> EdgeVec;
    EdgeVec.reserve(3 * NumTriangles());
    // This initial iteration will also fill the tri-edge adjacency
    m_T2E.resize(NumTriangles(), 3);
    for (int i = 0; i < NumTriangles(); ++i)
    {
        std::pair<int, int> e;
        std::pair<int, int> data;

        // Create the edge 01
        e.first = Triangles()(i, 0);
        e.second = Triangles()(i, 1);
        if (e.first > e.second)
            std::swap(e.first, e.second);
        data.first = i;
        data.second = EdgeVec.size();
        EdgeVec.emplace_back(e, data);
        // Triangle-edge adjacency
        m_T2E(i, 2) = data.second;

        // Create the edge 12
        e.first = Triangles()(i, 1);
        e.second = Triangles()(i, 2);
        if (e.first > e.second)
            std::swap(e.first, e.second);
        data.first = i;
        data.second = EdgeVec.size();
        EdgeVec.emplace_back(e, data);
        // Triangle-edge adjacency
        m_T2E(i, 0) = data.second;

        // Create the edge 20
        e.first = Triangles()(i, 2);
        e.second = Triangles()(i, 0);
        if (e.first > e.second)
            std::swap(e.first, e.second);
        data.first = i;
        data.second = EdgeVec.size();
        EdgeVec.emplace_back(e, data);
        // Triangle-edge adjacency
        m_T2E(i, 1) = data.second;
    }

    // Sort the edges to identify the uniques
    std::sort(EdgeVec.begin(), EdgeVec.end());
    int ELen = 0;
    std::pair<int, int> eCur = EdgeVec[0].first;
    for (int i = 1; i < EdgeVec.size(); ++i)
    {
        std::pair<int, int> e;
        std::pair<int, int> data;
        std::tie(e, data) = EdgeVec[i];
        if (e != eCur)
        {
            EdgeVec[++ELen].first = e;
            eCur = e;
        }
        for (int j = 0; j < 3; ++j)
        {
            if (m_T2E(data.first, j) == data.second)
            {
                m_T2E(data.first, j) = ELen;
                break;
            }
        }
    }
    ELen++;

    // Get the edges
    m_E.resize(ELen, 2);
    for (int i = 0; i < ELen; ++i)
    {
        m_E(i, 0) = EdgeVec[i].first.first;
        m_E(i, 1) = EdgeVec[i].first.second;
    }

    // Compute edge-triangle adjacency
    m_E2T.setConstant(NumEdges(), 2, -1);
    for (int i = 0; i < NumTriangles(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (m_E2T(m_T2E(i, j), 0) == -1)
                m_E2T(m_T2E(i, j), 0) = i;
            else
                m_E2T(m_T2E(i, j), 1) = i;
        }
    }

    // Compute triangle-triangle adjacency
    m_T2T.resize(NumTriangles(), 3);
    for (int i = 0; i < NumTriangles(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int t = m_E2T(m_T2E(i, j), 0);
            if (t == i)
                m_E2T(m_T2E(i, j), 1);
            m_T2T(i, j) = t;
        }
    }
}