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

#include <igl/cut_mesh.h>

#include <queue>


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

dfy::Mesh dfy::ManifoldMesh::CutEdges(const std::vector<int> &EdgeVec) const
{
    std::vector<bool> BFlags;
    BFlags.resize(NumEdges(), false);
    for (int e : EdgeVec)
        BFlags[e] = true;
    return CutEdges(BFlags);
}

dfy::Mesh dfy::ManifoldMesh::CutEdges(const Eigen::VectorXi &EdgeFlags) const
{
    std::vector<bool> BFlags;
    BFlags.reserve(NumEdges());
    for (int i = 0; i < NumEdges(); ++i)
        BFlags.emplace_back(EdgeFlags[i] > 0);
    return CutEdges(BFlags);
}

dfy::Mesh dfy::ManifoldMesh::CutEdges(const std::vector<bool> &EdgeFlags) const
{
    Eigen::MatrixXi TriFlags;
    TriFlags.setZero(NumTriangles(), 3);
    for (int i = 0; i < NumTriangles(); ++i)
    {
        for (int j = 0; j < 3; ++j)
            TriFlags(i, (j + 1) % 3) = EdgeFlags[TriEdgeAdj()(i, j)];
    }
    Eigen::MatrixXd VNew;
    Eigen::MatrixXi FNew;
    igl::cut_mesh(Vertices(), Triangles(), TriFlags, VNew, FNew);
    return dfy::Mesh(VNew, FNew);
}


void dfy::ManifoldMesh::InitManifoldMesh()
{
    // Create edges through a map
    std::map<std::pair<int, int>, int> EMap;
    // Use this to also fill tri-edge adjacency
    m_T2E.setConstant(NumTriangles(), 3, -1);
    for (int i = 0; i < NumTriangles(); ++i)
    {
        std::pair<int, int> e;
        for (int j = 0; j < 3; ++j)
        {
            e.first = Triangles()(i, j);
            e.second = Triangles()(i, (j + 1) % 3);
            if (e.second < e.first)
                std::swap(e.first, e.second);
            if (EMap.find(e) == EMap.end())
                EMap.emplace(e, (int)EMap.size());
            m_T2E(i, (j + 2) % 3) = EMap[e];
        }
    }

    // Get the edges
    m_E.resize(EMap.size(), 2);
    for (auto it : EMap)
    {
        m_E(it.second, 0) = it.first.first;
        m_E(it.second, 1) = it.first.second;
    }
    m_EL = (Vertices()(m_E(Eigen::all, 0), Eigen::all) - Vertices()(m_E(Eigen::all, 1), Eigen::all)).rowwise().norm();

    // Compute edge-triangle adjacency
    m_E2T.setConstant(NumEdges(), 2, -1);
    for (int i = 0; i < NumTriangles(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (m_T2E(i, j) == -1)
                continue;
            if (m_E2T(m_T2E(i, j), 0) == -1)
                m_E2T(m_T2E(i, j), 0) = i;
            else
                m_E2T(m_T2E(i, j), 1) = i;
        }
    }

    // Compute triangle-triangle adjacency
    m_T2T.setConstant(NumTriangles(), 3, -1);
    for (int i = 0; i < NumTriangles(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int e = m_T2E(i, j);
            if (e == -1)
                continue;
            int t = m_E2T(e, 0);
            if (t == i)
                t = m_E2T(e, 1);
            m_T2T(i, j) = t;
        }
    }
}


void dfy::ManifoldMesh::FlipTriangle(int i)
{
    dfy::Mesh::FlipTriangle(i);
    std::swap(m_T2E(i, 0), m_T2E(i, 1));
    std::swap(m_T2T(i, 0), m_T2T(i, 1));
}

void dfy::ManifoldMesh::FixFaceOrientation()
{
    std::vector<bool> Fixed;
    Fixed.resize(NumTriangles());

    std::queue<int> Q;
    Q.emplace(0);
    while (!Q.empty())
    {
        int CurFace = Q.front();
        Q.pop();

        if (Fixed[CurFace])
            continue;
        Fixed[CurFace] = true;

        Eigen::Vector3i Flip;
        Flip.setZero();
        for (int j = 0; j < 3; ++j)
        {
            int NextFace = TriTriAdj()(CurFace, j);
            if (NextFace < 0)
                continue;
            std::pair<int, int> CurE{ Triangles()(CurFace, (j + 1) % 3), 
                                      Triangles()(CurFace, (j + 2) % 3) };
            std::pair<int, int> NextE;
            for (int k = 0; k < 3; ++k)
            {
                if (TriTriAdj()(NextFace, k) != CurFace)
                    continue;
                NextE.first = Triangles()(NextFace, (k + 1) % 3);
                NextE.second = Triangles()(NextFace, (k + 2) % 3);
            }
            // If shared edge has same orientation, we have to flip the triangle
            if (CurE == NextE)
                Flip[j] = 1;
            Q.emplace(NextFace);
        }

        for (int j = 0; j < 3; ++j)
        {
            if (Flip[j] == 0)
                continue;
            FlipTriangle(TriTriAdj()(CurFace, j));
        }
    }

    // Try estimating inside-outside to check if all faces are oriented outside
    Eigen::RowVector3d Center = FaceBarycs().colwise().mean();
    Eigen::MatrixXd GoInside = FaceBarycs() - Center.replicate(NumTriangles(), 1);
    Eigen::VectorXi DotProds(NumTriangles());
    for (int i = 0; i < NumTriangles(); ++i)
        DotProds(i) = GoInside.row(i).dot(Normals().row(i)) < 0 ? 1 : -1;
    if (DotProds.sum() < 0)
        return;
    for (int i = 0; i < NumTriangles(); ++i)
        FlipTriangle(i);
}