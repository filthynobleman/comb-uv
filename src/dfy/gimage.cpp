/**
 * @file        gimage.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-14
 */
#include <dfy/gimage.hpp>

#include <dfy/utils.hpp>

#include <igl/barycentric_coordinates.h>

dfy::GImage::GImage(const dfy::Mesh &M, 
                    int Width, 
                    int Height)
    : m_Mesh(M), m_Width(Width), m_Height(Height)
{
    m_Data = (Eigen::Vector3d*)std::calloc(Width * Height, sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot create GImage.");
}

dfy::GImage::GImage(const dfy::Mesh &M, int Size)
    : dfy::GImage(M, Size, Size) { }

dfy::GImage::GImage(const dfy::GImage &GI)
    : m_Mesh(GI.m_Mesh)
{
    m_Width = GI.m_Width;
    m_Height = GI.m_Height;
    size_t BufSize = m_Width * m_Height * sizeof(Eigen::Vector3d);
    m_Data = (Eigen::Vector3d*)std::malloc(BufSize);
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot create GImage.");
    m_Data = (Eigen::Vector3d*)std::memcpy(m_Data, GI.m_Data, BufSize);
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot create GImage.");
}

dfy::GImage::GImage(dfy::GImage &&GI)
    : m_Mesh(GI.m_Mesh)
{
    m_Width = GI.m_Width;
    m_Height = GI.m_Height;
    m_Data = GI.m_Data;
    GI.m_Data = nullptr;
}

dfy::GImage::~GImage()
{
    delete m_Data;
}

int dfy::GImage::GetWidth() const { return m_Width; }
int dfy::GImage::GetHeight() const { return m_Height; }

const Eigen::Vector3d &dfy::GImage::GetPixel(int i, int j) const
{
    return m_Data[j * GetWidth() + i];
}

void dfy::GImage::GetPixel(int i, int j, 
                           Eigen::Vector3d &rgb) const
{
    rgb = GetPixel(i, j);
}

void dfy::GImage::GetPixel(int i, int j, 
                           unsigned char &r, unsigned char &g, unsigned char &b) const
{
    const Eigen::Vector3d& rgb = GetPixel(i, j);
    r = std::round(std::min(1.0, std::max(0.0, rgb.x())) * 255.0);
    g = std::round(std::min(1.0, std::max(0.0, rgb.y())) * 255.0);
    b = std::round(std::min(1.0, std::max(0.0, rgb.z())) * 255.0);
}

#include <iostream>
#include <random>

void dfy::GImage::Compute(const Eigen::MatrixXd &UV, 
                          const Eigen::MatrixXi &TUV)
{
    // Normalize mesh coordinates
    Eigen::MatrixXd RGB = m_Mesh.Vertices();
    for (int i = 0; i < 3; ++i)
        RGB.col(i) = RGB.col(i).array() - RGB.col(i).minCoeff();
    RGB /= RGB.maxCoeff();
    for (int t = 0; t < TUV.rows(); ++t)
    {
        // Get triangle
        Eigen::Matrix<double, 3, 2> Tri = UV(TUV(t, Eigen::all).transpose(), Eigen::all);
        // Get bounding box of triangle
        Eigen::Vector2d dBBL = Tri.colwise().minCoeff();
        Eigen::Vector2d dBTR = Tri.colwise().maxCoeff();
        Eigen::Vector2i BBL, BTR;
        BBL[0] = std::floor(dBBL[0] * GetWidth());
        BBL[1] = std::floor(dBBL[1] * GetHeight());
        BTR[0] = std::ceil(dBTR[0] * GetWidth());
        BTR[1] = std::ceil(dBTR[1] * GetHeight());
        // Get triangle's coordinate in pixel space
        // Tri.col(0) = Tri.col(0).array() * GetWidth();
        // Tri.col(1) = Tri.col(1).array() * GetHeight();
        
        // Iterate over the pixels inside the bounding box
        for (int j = BBL.y(); j < BTR.y(); ++j)
        {
            for (int i = BBL.x(); i < BTR.x(); ++i)
            {
                // Get the pixel coordinates
                // Eigen::Vector2d ij{ i, j };
                Eigen::RowVector2d ij{ (i + 1e-9) / (GetWidth() - (1 - 2e-9)), (j + 1e-9) / (GetHeight() - (1 - 2e-9)) };
                Eigen::RowVector3d Lambda;
                igl::barycentric_coordinates(ij, Tri.row(0), Tri.row(1), Tri.row(2), Lambda);
                if ((Lambda.array() > -1e-6).all())
                {
                    Eigen::Matrix3d T3D = RGB(m_Mesh.Triangles()(t, Eigen::all), Eigen::all);
                    m_Data[j * GetWidth() + i] = T3D.transpose() * Lambda.transpose();
                }
                // // Get barycentric coordinates of pixel
                // Eigen::Matrix2d T;
                // T.row(0) = Tri.row(0) - Tri.row(2);
                // T.row(1) = Tri.row(1) - Tri.row(2);
                // ij = (T.inverse() * (ij - Tri.row(2).transpose())).eval();
                // // If outside, ignore
                // // if (ij[0] < 0 || ij[1] < 0 || ij.sum() > 1)
                // //     continue;
                // // m_Data[j * GetWidth() + i].setOnes();
                // // Get the 3D coordinates using barycentric interpolation inside the triangle
                // m_Data[j * GetWidth() + i] = ij[0] * m_Mesh.Vertices()(m_Mesh.Triangles()(i, 0), Eigen::all);
                // m_Data[j * GetWidth() + i] += ij[1] * m_Mesh.Vertices()(m_Mesh.Triangles()(i, 1), Eigen::all);
                // m_Data[j * GetWidth() + i] += (1.0 - ij.sum()) * m_Mesh.Vertices()(m_Mesh.Triangles()(i, 2), Eigen::all);
            }
        }
    }
}