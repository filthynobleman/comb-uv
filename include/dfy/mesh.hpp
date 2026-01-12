/**
 * @file        mesh.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-09
 */
#pragma once

#include <vector>
#include <string>

#include <Eigen/Dense>


namespace dfy
{

enum SubMeshAccessType {
    BY_VERTICES = 1,
    BY_VERTICES_ALL = BY_VERTICES | (BY_VERTICES << 1),
    BY_VERTICES_ANY = BY_VERTICES,
    BY_TRIANGLES = BY_VERTICES << 2
};

enum NormalType {
    AT_VERTICES,
    AT_TRIANGLES
};
    
class Mesh
{
private:
    // Mesh vertices represented as a #V-by-3 double matrix
    Eigen::MatrixXd m_V;

    // Mesh triangles represented as a #F-by-3 integer matrix
    Eigen::MatrixXi m_F;

    // Mesh triangle areas
    Eigen::VectorXd m_Areas;

    // Mesh triangle angles
    Eigen::MatrixXd m_Angles;

    // Mesh triangle barycenters
    Eigen::MatrixXd m_Barycs;

    // Mesh normals at vertices
    Eigen::MatrixXd m_NV;

    // Mesh normals at triangles
    Eigen::MatrixXd m_NF;

    // Initialization procedure
    void InitMesh();

public:
    Mesh(const Eigen::MatrixXd& Verts,
         const Eigen::MatrixXi& Tris);
    Mesh(const std::string& Filename);
    Mesh(const dfy::Mesh& M);
    Mesh(dfy::Mesh&& M);
    virtual dfy::Mesh operator=(const dfy::Mesh& M);
    virtual dfy::Mesh operator=(dfy::Mesh&& M);
    ~Mesh();

    int NumVertices() const;
    int NumTriangles() const;

    const Eigen::MatrixXd& Vertices() const;
    const Eigen::MatrixXi& Triangles() const;
    const Eigen::MatrixXd& Normals(dfy::NormalType Type = dfy::NormalType::AT_TRIANGLES) const;
    const Eigen::MatrixXd& VertNormals() const;
    const Eigen::MatrixXd& FaceNormals() const;
    const Eigen::VectorXd& FaceAreas() const;
    const Eigen::MatrixXd& FaceAngles() const;
    const Eigen::MatrixXd& FaceBarycs() const;


    dfy::Mesh SubMesh(const std::vector<int>& Indices,
                      dfy::SubMeshAccessType AccType,
                      std::vector<int>& VIdx,
                      std::vector<int>& TIdx) const;
    dfy::Mesh SubMesh(const Eigen::VectorXi& Indices,
                      dfy::SubMeshAccessType AccType,
                      std::vector<int>& VIdx,
                      std::vector<int>& TIdx) const;
    dfy::Mesh SubMesh(const int* const Indices,
                      int nIndices,
                      dfy::SubMeshAccessType AccType,
                      std::vector<int>& VIdx,
                      std::vector<int>& TIdx) const;

    dfy::Mesh SubMesh(const std::vector<int>& Indices,
                      dfy::SubMeshAccessType AccType) const;
    dfy::Mesh SubMesh(const Eigen::VectorXi& Indices,
                      dfy::SubMeshAccessType AccType) const;
    dfy::Mesh SubMesh(const int* const Indices,
                      int nIndices,
                      dfy::SubMeshAccessType AccType) const;
};

} // namespace dfy
