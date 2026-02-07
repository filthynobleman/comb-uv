/**
 * @file        mesh.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-09
 */
#include <dfy/mesh.hpp>
#include <dfy/utils.hpp>

#include <igl/read_triangle_mesh.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/barycenter.h>



dfy::Mesh::Mesh(const Eigen::MatrixXd &Verts, 
                const Eigen::MatrixXi &Tris)
{
    m_V = Verts;
    m_F = Tris;
    InitMesh();
}

dfy::Mesh::Mesh(const std::string &Filename)
{
    igl::read_triangle_mesh(Filename, m_V, m_F);
    InitMesh();
}

dfy::Mesh::Mesh(const dfy::Mesh &M)
{
    m_V = M.m_V;
    m_F = M.m_F;
    m_Areas = M.m_Areas;
    m_Angles = M.m_Angles;
    m_Barycs = M.m_Barycs;
    m_NV = M.m_NV;
    m_NF = M.m_NF;
}

dfy::Mesh::Mesh(dfy::Mesh &&M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_Areas = std::move(M.m_Areas);
    m_Angles = std::move(M.m_Angles);
    m_Barycs = std::move(M.m_Barycs);
    m_NV = std::move(M.m_NV);
    m_NF = std::move(M.m_NF);
}

dfy::Mesh dfy::Mesh::operator=(const dfy::Mesh &M)
{
    m_V = M.m_V;
    m_F = M.m_F;
    m_Areas = M.m_Areas;
    m_Angles = M.m_Angles;
    m_Barycs = M.m_Barycs;
    m_NV = M.m_NV;
    m_NF = M.m_NF;
    return *this;
}

dfy::Mesh dfy::Mesh::operator=(dfy::Mesh &&M)
{
    m_V = std::move(M.m_V);
    m_F = std::move(M.m_F);
    m_Areas = std::move(M.m_Areas);
    m_Angles = std::move(M.m_Angles);
    m_Barycs = std::move(M.m_Barycs);
    m_NV = std::move(M.m_NV);
    m_NF = std::move(M.m_NF);
    return *this;
}

dfy::Mesh::~Mesh() { }

void dfy::Mesh::InitMesh()
{
    igl::per_vertex_normals(m_V, m_F, m_NV);
    igl::per_face_normals(m_V, m_F, m_NF);
    igl::doublearea(m_V, m_F, m_Areas);
    m_Areas = 0.5 * m_Areas.array();
    igl::internal_angles(m_V, m_F, m_Angles);
    igl::barycenter(m_V, m_F, m_Barycs);
}

int dfy::Mesh::NumVertices() const { return m_V.rows(); }
int dfy::Mesh::NumTriangles() const { return m_F.rows(); }

void dfy::Mesh::FlipTriangle(int i)
{
    std::swap(m_F(i, 0), m_F(i, 1));
    std::swap(m_Angles(i, 0), m_Angles(i, 1));
    m_NF.row(i) *= -1;
}


const Eigen::MatrixXd& dfy::Mesh::Vertices() const { return m_V; }
const Eigen::MatrixXi& dfy::Mesh::Triangles() const { return m_F; }
const Eigen::MatrixXd &dfy::Mesh::VertNormals() const { return m_NV; }
const Eigen::MatrixXd &dfy::Mesh::FaceNormals() const { return m_NF; }

const Eigen::MatrixXd &dfy::Mesh::Normals(dfy::NormalType Type) const
{
    switch (Type)
    {
    case dfy::NormalType::AT_VERTICES:
        return m_NV;
    
    case dfy::NormalType::AT_TRIANGLES:
        return m_NF;

    default:
        throw std::runtime_error("Given unsupported normal type.");
    }
}

const Eigen::VectorXd &dfy::Mesh::FaceAreas() const { return m_Areas; }
const Eigen::MatrixXd &dfy::Mesh::FaceAngles() const { return m_Angles; }
const Eigen::MatrixXd &dfy::Mesh::FaceBarycs() const { return m_Barycs; }

dfy::Mesh dfy::Mesh::SubMesh(const int *const Indices, int nIndices, 
                             dfy::SubMeshAccessType AccType, 
                             std::vector<int> &VIdx, 
                             std::vector<int> &TIdx) const
{
    Eigen::VectorXi Potential;
    Eigen::VectorXi TriPotential;
    // If access type is by vertices, we must obtain the triangles
    if ((AccType & dfy::SubMeshAccessType::BY_VERTICES) != 0)
    {
        Potential.setZero(NumVertices());
        for (int k = 0; k < nIndices; ++k)
            Potential[Indices[k]] = 1;
        TriPotential.setZero(NumTriangles());
        for (int i = 0; i < NumTriangles(); ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (Potential(Triangles()(i, j)) == 0)
                    continue;
                TriPotential[i]++;
            }
        }
        int Thresh = 1;
        if (AccType == dfy::SubMeshAccessType::BY_VERTICES_ALL)
            Thresh = 3;
        for (int i = 0; i < NumTriangles(); ++i)
            TriPotential[i] = TriPotential[i] >= Thresh ? 1 : 0;
    }
    // Otherwise, we already have the triangles
    else
    {
        TriPotential.setZero(NumTriangles());
        for (int k = 0; k < nIndices; ++k)
            TriPotential[Indices[k]] = 1;
        Potential.setZero(NumVertices());
    }

    // Get the right vertex potential
    for (int i = 0; i < NumTriangles(); ++i)
    {
        if (TriPotential[i] == 0)
            continue;
        for (int j = 0; j < 3; ++j)
            Potential[Triangles()(i, j)] = 1;
    }

    // Extract vertex indices
    VIdx.clear();
    VIdx.reserve(Potential.sum());
    std::map<int, int> VMap;
    for (int i = 0; i < NumVertices(); ++i)
    {
        if (Potential[i] == 0)
            continue;
        
        VMap.emplace(i, (int)VIdx.size());
        VIdx.emplace_back(i);
    }

    // Make triangles
    TIdx.clear();
    TIdx.reserve(TriPotential.sum());
    for (int i = 0; i < NumTriangles(); ++i)
    {
        if (TriPotential[i] == 0)
            continue;
        
        TIdx.emplace_back(i);
    }

    Eigen::MatrixXi Tris = m_F(TIdx, Eigen::all);
    for (int i = 0; i < Tris.rows(); ++i)
    {
        for (int j = 0; j < 3; ++j)
            Tris(i, j) = VMap[Tris(i, j)];
    }

    // Return submesh
    return dfy::Mesh(m_V(VIdx, Eigen::all), Tris);
}

dfy::Mesh dfy::Mesh::SubMesh(const int *const Indices, int nIndices, 
                             dfy::SubMeshAccessType AccType) const
{
    std::vector<int> VIdx;
    std::vector<int> TIdx;
    return SubMesh(Indices, nIndices, AccType, VIdx, TIdx);
}

dfy::Mesh dfy::Mesh::SubMesh(const std::vector<int> &Indices, 
                             dfy::SubMeshAccessType AccType, 
                             std::vector<int> &VIdx, 
                             std::vector<int> &TIdx) const
{
    return SubMesh(Indices.data(), Indices.size(), AccType, VIdx, TIdx);
}

dfy::Mesh dfy::Mesh::SubMesh(const std::vector<int> &Indices, 
                             dfy::SubMeshAccessType AccType) const
{
    return SubMesh(Indices.data(), Indices.size(), AccType);
}

dfy::Mesh dfy::Mesh::SubMesh(const Eigen::VectorXi &Indices, 
                             dfy::SubMeshAccessType AccType, 
                             std::vector<int> &VIdx, 
                             std::vector<int> &TIdx) const
{
    return SubMesh(Indices.data(), Indices.rows(), AccType, VIdx, TIdx);
}

dfy::Mesh dfy::Mesh::SubMesh(const Eigen::VectorXi &Indices, 
                             dfy::SubMeshAccessType AccType) const
{
    return SubMesh(Indices.data(), Indices.rows(), AccType);
}