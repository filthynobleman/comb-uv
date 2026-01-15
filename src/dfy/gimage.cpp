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
    // Check each pixel has been set
    std::vector<bool> PixSet;
    PixSet.resize(GetWidth() * GetHeight(), false);
    for (int t = 0; t < TUV.rows(); ++t)
    {
        // Get triangle
        Eigen::Matrix<double, 3, 2> Tri = UV(TUV(t, Eigen::all).transpose(), Eigen::all);
        // Get bounding box of triangle
        Eigen::Vector2d dBBL = Tri.colwise().minCoeff().cwiseMax(0.0);
        Eigen::Vector2d dBTR = Tri.colwise().maxCoeff().cwiseMin(1.0);
        Eigen::Vector2i BBL, BTR;
        BBL[0] = std::floor(dBBL[0] * GetWidth());
        BBL[1] = std::floor(dBBL[1] * GetHeight());
        BTR[0] = std::ceil(dBTR[0] * GetWidth());
        BTR[1] = std::ceil(dBTR[1] * GetHeight());
        
        // Iterate over the pixels inside the bounding box
        for (int j = BBL.y(); j < BTR.y(); ++j)
        {
            for (int i = BBL.x(); i < BTR.x(); ++i)
            {
                // Get the pixel coordinates
                Eigen::RowVector2d ij{ (i + 1e-6) / (GetWidth() - (1 - 2e-6)), (j + 1e-6) / (GetHeight() - (1 - 2e-6)) };
                Eigen::RowVector3d Lambda;
                igl::barycentric_coordinates(ij, Tri.row(0), Tri.row(1), Tri.row(2), Lambda);
                if ((Lambda.array() > -1e-6).all())
                {
                    Eigen::Matrix3d T3D = RGB(m_Mesh.Triangles()(t, Eigen::all), Eigen::all);
                    m_Data[j * GetWidth() + i] = Lambda * T3D;
                    PixSet[j * GetWidth() + i] = true;
                }
            }
        }
    }

    // Fill unset pixels with closest triangle
    for (int j = 0; j < GetHeight(); ++j)
    {
        for (int i = 0; i < GetWidth(); ++i)
        {
            if (PixSet[j * GetWidth() + i])
                continue;

            // Get the pixel coordinates
            Eigen::RowVector2d ij{ (i + 1e-6) / (GetWidth() - (1 - 2e-6)), 
                                   (j + 1e-6) / (GetHeight() - (1 - 2e-6)) };
            
            // Find closest triangle
            int Closest = -1;
            double Error = std::numeric_limits<double>::infinity();
            for (int t = 0; t < TUV.rows(); ++t)
            {
                // Get triangle
                Eigen::Matrix<double, 3, 2> Tri = UV(TUV(t, Eigen::all).transpose(), Eigen::all);
                // Get barycentric coords
                Eigen::RowVector3d Lambda;
                igl::barycentric_coordinates(ij, Tri.row(0), Tri.row(1), Tri.row(2), Lambda);
                Lambda = Lambda.cwiseMax(0);
                Lambda /= Lambda.sum();
                double err = (ij - Lambda * Tri).eval().squaredNorm();
                if (err < Error)
                {
                    Error = err;
                    Closest = t;
                }
            }

            // Apply closest triangle
            int t = Closest;
            // Get triangle
            Eigen::Matrix<double, 3, 2> Tri = UV(TUV(t, Eigen::all).transpose(), Eigen::all);
            // Get barycentric coords
            Eigen::RowVector3d Lambda;
            igl::barycentric_coordinates(ij, Tri.row(0), Tri.row(1), Tri.row(2), Lambda);
            Lambda = Lambda.cwiseMax(0);
            Lambda /= Lambda.sum();
            Eigen::Matrix3d T3D = RGB(m_Mesh.Triangles()(t, Eigen::all), Eigen::all);
            m_Data[j * GetWidth() + i] = Lambda * T3D;
        }
    }
}