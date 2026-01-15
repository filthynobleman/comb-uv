/**
 * @file        quadmesh.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-15
 */
#include <dfy/quadmesh.hpp>

dfy::QuadMesh::QuadMesh(const Eigen::MatrixXd &Verts, 
                        const Eigen::MatrixXi &Quads)
{
    m_V = Verts;
    m_F = Quads;
}

dfy::QuadMesh::QuadMesh(const Eigen::MatrixXd &Grid, 
                        int Width, int Height)
{
    m_V = Grid;
    m_F.resize((Width - 1) * (Height - 1), 4);
    for (int j = 0; j < Height - 1; ++j)
    {
        for (int i = 0; i < Width - 1; ++i)
        {
            int V00 = j * Width + i;
            int V01 = V00 + 1;
            int V10 = V00 + Width;
            int V11 = V10 + 1;
            int Row = j * (Width - 1) + i;
            m_F.row(Row) = Eigen::Vector4i{ V00, V01, V11, V10 };
        }
    }
}

dfy::QuadMesh::QuadMesh(const dfy::QuadMesh &M)
{
    m_V = M.m_V;
    m_F = M.m_F;
}

dfy::QuadMesh::QuadMesh(dfy::QuadMesh &&M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
}

dfy::QuadMesh& dfy::QuadMesh::operator=(const dfy::QuadMesh &M)
{
    m_V = M.m_V;
    m_F = M.m_F;
    return *this;
}

dfy::QuadMesh& dfy::QuadMesh::operator=(dfy::QuadMesh &&M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    return *this;
}

dfy::QuadMesh::~QuadMesh() { }

int dfy::QuadMesh::NumVertices() const { return m_V.rows(); }
int dfy::QuadMesh::NumQuads() const { return m_F.rows(); }
const Eigen::MatrixXd &dfy::QuadMesh::Vertices() const { return m_V; }
const Eigen::MatrixXi &dfy::QuadMesh::Quads() const { return m_F; }