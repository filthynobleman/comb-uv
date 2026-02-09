/**
 * @file        manifoldmesh.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#pragma once

#include <dfy/mesh.hpp>


namespace dfy
{
    
class ManifoldMesh : public dfy::Mesh
{
private:
    // Mesh edges
    Eigen::MatrixXi m_E;
    
    // Mesh edge lengths
    Eigen::VectorXd m_EL;

    // Edge-to-triangle adjacency
    Eigen::MatrixXi m_E2T;

    // Triangle-to-edge adjqcency
    Eigen::MatrixXi m_T2E;

    // Triangle-to-triangle adjacency
    Eigen::MatrixXi m_T2T;

    // Inizialization procedure
    void InitManifoldMesh();


public:
    ManifoldMesh(const Eigen::MatrixXd& Verts,
                 const Eigen::MatrixXi& Tris);
    ManifoldMesh(const std::string& Filename);
    ManifoldMesh(const dfy::ManifoldMesh& M);
    ManifoldMesh(dfy::ManifoldMesh&& M);
    virtual dfy::ManifoldMesh& operator=(const dfy::ManifoldMesh& M);
    virtual dfy::ManifoldMesh& operator=(dfy::ManifoldMesh&& M);
    ~ManifoldMesh();

    int NumEdges() const;
    const Eigen::MatrixXi& Edges() const;
    const Eigen::VectorXd& EdgeLengths() const;
    Eigen::MatrixXd EdgeCenters() const;
    const Eigen::MatrixXi& EdgeTriAdj() const;
    const Eigen::MatrixXi& TriEdgeAdj() const;
    const Eigen::MatrixXi& TriTriAdj() const;

    dfy::Mesh CutEdges(const std::vector<int>& EdgeVec) const;
    dfy::Mesh CutEdges(const std::vector<bool>& EdgeFlags) const;
    dfy::Mesh CutEdges(const Eigen::VectorXi& EdgeFlags) const;

    virtual void FlipTriangle(int i) override;
    void FixFaceOrientation();

    bool IsDisk() const;
};


} // namespace dfy
