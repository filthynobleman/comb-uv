/**
 * @file        gimage.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-14
 */
#pragma once

#include <dfy/mesh.hpp>
#include <dfy/quadmesh.hpp>



namespace dfy
{
    
class GImage
{
private:
    const dfy::Mesh& m_Mesh;

    Eigen::Vector3d* m_Data;
    int m_Width;
    int m_Height;

public:
    GImage(const dfy::Mesh& M,
           int Width,
           int Height);
    GImage(const dfy::Mesh& M,
           int Size);
    GImage(const dfy::GImage& GI);
    GImage(dfy::GImage&& GI);
    ~GImage();

    int GetWidth() const;
    int GetHeight() const;

    const Eigen::Vector3d& GetPixel(int i, int j) const;
    void GetPixel(int i, int j,
                  Eigen::Vector3d& rgb) const;
    void GetPixel(int i, int j,
                  unsigned char& r, unsigned char& g, unsigned char& b) const;


    void Compute(const Eigen::MatrixXd& UV,
                 const Eigen::MatrixXi& TUV);

    dfy::QuadMesh AsQuadMesh(const dfy::Mesh& M, int URes, int VRes, bool Weld = false) const;
    dfy::QuadMesh AsQuadMesh(const dfy::Mesh& M, int Res, bool Weld = false) const;
};


} // namespace dfy
