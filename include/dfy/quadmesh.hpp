/**
 * @file        quadmesh.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-15
 */
#pragma once


#include <dfy/utils.hpp>


namespace dfy
{
    
class QuadMesh
{
private:
    Eigen::MatrixXd m_V;
    Eigen::MatrixXi m_F;

public:
    QuadMesh(const Eigen::MatrixXd& Verts,
             const Eigen::MatrixXi& Quads);
    QuadMesh(const Eigen::MatrixXd& Grid,
             int Width, int Height);
    QuadMesh(const dfy::QuadMesh& M);
    QuadMesh(dfy::QuadMesh&& M);
    dfy::QuadMesh& operator=(const dfy::QuadMesh& M);
    dfy::QuadMesh& operator=(dfy::QuadMesh&& M);
    ~QuadMesh();

    int NumVertices() const;
    int NumQuads() const;

    const Eigen::MatrixXd& Vertices() const;
    const Eigen::MatrixXi& Quads() const;
};

} // namespace dfy
