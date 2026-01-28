/**
 * @file        imgsampler.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-15
 */
#include <dfy/imgsampler.hpp>

#include <stb_image.h>


dfy::ImageSampler::ImageSampler(const std::string &Filename)
{
    stbi_set_flip_vertically_on_load(true);
    int nCh;
    unsigned char* img = stbi_load(Filename.c_str(), &m_Width, &m_Height, &nCh, 3);
    if (img == nullptr)
        throw std::runtime_error("Cannot load image " + Filename);
    

    size_t BufSize = m_Width * m_Height;
    m_Data = (Eigen::Vector3d*)std::malloc(BufSize * sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot allocate image sampler.");
    
    for (int j = 0; j < m_Height; ++j)
    {
        for (int i = 0; i < m_Width; ++i)
        {
            for (int k = 0; k < 3; ++k)
            {
                size_t Idx = 3 * (m_Width * j + i) + k;
                double Val = img[Idx] / 255.0;
                m_Data[m_Width * j + i][k] = Val;
            }
        }
    }
}

dfy::ImageSampler::ImageSampler(const dfy::GImage &Img)
{
    m_Width = Img.GetWidth();
    m_Height = Img.GetHeight();
    size_t BufSize = m_Width * m_Height;
    m_Data = (Eigen::Vector3d*)std::malloc(BufSize * sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot allocate image sampler.");

    for (int j = 0; j < m_Height; ++j)
    {
        for (int i = 0; i < m_Width; ++i)
            Img.GetPixel(i, j, m_Data[j * m_Width + i]);
    }
}

dfy::ImageSampler::ImageSampler(const dfy::ImageSampler &IS)
{
    m_Width = IS.m_Width;
    m_Height = IS.m_Height;
    size_t BufSize = m_Width * m_Height;
    m_Data = (Eigen::Vector3d*)std::malloc(BufSize * sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot allocate image sampler.");
    m_Data = (Eigen::Vector3d*)std::memcpy(m_Data, IS.m_Data, BufSize * sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot allocate image sampler.");
}

dfy::ImageSampler::ImageSampler(dfy::ImageSampler &&IS)
{
    m_Width = IS.m_Width;
    m_Height = IS.m_Height;
    m_Data = IS.m_Data;
    IS.m_Data = nullptr;
}

dfy::ImageSampler& dfy::ImageSampler::operator=(const dfy::ImageSampler &IS)
{
    m_Width = IS.m_Width;
    m_Height = IS.m_Height;
    size_t BufSize = m_Width * m_Height;
    m_Data = (Eigen::Vector3d*)std::malloc(BufSize * sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot allocate image sampler.");
    m_Data = (Eigen::Vector3d*)std::memcpy(m_Data, IS.m_Data, BufSize * sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot allocate image sampler.");

    return *this;
}

dfy::ImageSampler& dfy::ImageSampler::operator=(dfy::ImageSampler &&IS)
{
    m_Width = IS.m_Width;
    m_Height = IS.m_Height;
    m_Data = IS.m_Data;
    IS.m_Data = nullptr;

    return *this;
}

dfy::ImageSampler::~ImageSampler()
{
    delete m_Data;
}

int dfy::ImageSampler::GetWidth() const { return m_Width; }
int dfy::ImageSampler::GetHeight() const { return m_Height; }

const Eigen::Vector3d &dfy::ImageSampler::GetPixel(int i, int j) const
{
    return m_Data[j * m_Width + i];
}

Eigen::Vector3d dfy::ImageSampler::Sample(const Eigen::Vector2d &uv,
                                          dfy::ImgInterpType InterpType) const
{
    return Sample(uv[0], uv[1], InterpType);
}

double Q1(double x) { return 0.5 * (x * (x * (2 - x) - 1)); }
double Q2(double x) { return 0.5 * (2 + x * x * (3 * x - 5)); }
double Q3(double x) { return 0.5 * (x * (1 + x * (4 - 3 * x))); }
double Q4(double x) { return 0.5 * (x * x * (x - 1)); }

double dQ1(double x) { return 0.5 * (4*x - 3*x*x - 1); }
double dQ2(double x) { return 0.5 * (9*x*x - 10*x); }
double dQ3(double x) { return 0.5 * (1 + 8*x - 9*x*x); }
double dQ4(double x) { return 0.5 * (3*x*x - 2*x); }

Eigen::Vector3d dfy::ImageSampler::Sample(double u, double v, 
                                          dfy::ImgInterpType InterpType) const
{
    // Clamp uv to [0, 1]
    u = std::min(1.0, std::max(0.0, u));
    v = std::min(1.0, std::max(0.0, v));
    // Get surrounding pixels
    u *= GetWidth() - 1;
    int u1 = std::floor(u);
    int u2 = std::ceil(u);
    u -= u1;
    v *= GetHeight() - 1;
    int v1 = std::floor(v);
    int v2 = std::ceil(v);
    v -= v1;
    if (InterpType == dfy::ImgInterpType::LINEAR)
    {
        Eigen::Vector3d D1 = u * GetPixel(u1, v1) + (1 - u) * GetPixel(u2, v1);
        Eigen::Vector3d D2 = u * GetPixel(u1, v2) + (1 - u) * GetPixel(u2, v2);
        return v * D2 + (1 - v) * D1;
    }
    else if (InterpType == dfy::ImgInterpType::CUBIC)
    {
        int u0 = std::max(0, u1 - 1);
        int u3 = std::min(GetWidth() - 1, u2 + 1);
        int v0 = std::max(0, v1 - 1);
        int v3 = std::min(GetHeight() - 1, v2 + 1);
        double uq0 = Q1(u);
        double uq1 = Q2(u);
        double uq2 = Q3(u);
        double uq3 = Q4(u);
        Eigen::Vector3d D0 = uq0 * GetPixel(u0, v0) +
                             uq1 * GetPixel(u1, v0) +
                             uq2 * GetPixel(u2, v0) +
                             uq3 * GetPixel(u3, v0);
        Eigen::Vector3d D1 = uq0 * GetPixel(u0, v1) +
                             uq1 * GetPixel(u1, v1) +
                             uq2 * GetPixel(u2, v1) +
                             uq3 * GetPixel(u3, v1);
        Eigen::Vector3d D2 = uq0 * GetPixel(u0, v2) +
                             uq1 * GetPixel(u1, v2) +
                             uq2 * GetPixel(u2, v2) +
                             uq3 * GetPixel(u3, v2);
        Eigen::Vector3d D3 = uq0 * GetPixel(u0, v3) +
                             uq1 * GetPixel(u1, v3) +
                             uq2 * GetPixel(u2, v3) +
                             uq3 * GetPixel(u3, v3);
        return Q1(v) * D0 + Q2(v) * D1 + Q3(v) * D2 + Q4(v) * D3;
    }
    else
        throw std::runtime_error("Undefined interpolation scheme.");
}

Eigen::Matrix<double, 3, 2> dfy::ImageSampler::Jacobian(const Eigen::Vector2d &uv, 
                                                        dfy::ImgInterpType InterpType) const
{
    return Jacobian(uv.x(), uv.y(), InterpType);
}

Eigen::Matrix<double, 3, 2> dfy::ImageSampler::Jacobian(double u, double v, 
                                                        dfy::ImgInterpType InterpType) const
{
    Eigen::Matrix<double, 3, 2> J;
    const double h = 1e-6;
    J.col(0) = (Sample(u + h, v, InterpType) - Sample(u - h, v, InterpType)) / (2 * h);
    J.col(1) = (Sample(u, v + h, InterpType) - Sample(u, v - h, InterpType)) / (2 * h);
    return J;

    // // Clamp uv to [0, 1]
    // u = std::min(1.0, std::max(0.0, u));
    // v = std::min(1.0, std::max(0.0, v));
    // // Get surrounding pixels
    // u *= GetWidth() - 1;
    // int u1 = std::floor(u);
    // int u2 = std::ceil(u);
    // u -= u1;
    // v *= GetHeight() - 1;
    // int v1 = std::floor(v);
    // int v2 = std::ceil(v);
    // v -= v1;
    // if (InterpType == dfy::ImgInterpType::LINEAR)
    // {
    //     Eigen::Matrix<double, 3, 2> J;
    //     Eigen::Vector3d D1;
    //     Eigen::Vector3d D2;
    //     // Get horizontal extremes at v
    //     D1 = v * GetPixel(u1, v1) + (1 - v) * GetPixel(u1, v2);
    //     D2 = v * GetPixel(u2, v1) + (1 - v) * GetPixel(u2, v2);
    //     J.col(0) = D2 - D1;
    //     // Get horizontal extremes at u
    //     D1 = u * GetPixel(u1, v1) + (1 - u) * GetPixel(u2, v1);
    //     D2 = u * GetPixel(u1, v2) + (1 - u) * GetPixel(u2, v2);
    //     J.col(1) = D2 - D1;
    //     return J;
    // }
    // else if (InterpType == dfy::ImgInterpType::CUBIC)
    // {
    //     Eigen::Matrix<double, 3, 2> J;
    //     Eigen::Vector3d D0;
    //     Eigen::Vector3d D1;
    //     Eigen::Vector3d D2;
    //     Eigen::Vector3d D3;
    //     int u0 = std::max(0, u1 - 1);
    //     int u3 = std::min(GetWidth() - 1, u2 + 1);
    //     int v0 = std::max(0, v1 - 1);
    //     int v3 = std::min(GetHeight() - 1, v2 + 1);

    //     double vq0 = Q1(v);
    //     double vq1 = Q2(v);
    //     double vq2 = Q3(v);
    //     double vq3 = Q4(v);
    //     D0 = vq0 * GetPixel(u0, v0) +
    //          vq1 * GetPixel(u0, v1) +
    //          vq2 * GetPixel(u0, v2) +
    //          vq3 * GetPixel(u0, v3);
    //     D1 = vq0 * GetPixel(u1, v0) +
    //          vq1 * GetPixel(u1, v1) +
    //          vq2 * GetPixel(u1, v2) +
    //          vq3 * GetPixel(u1, v3);
    //     D2 = vq0 * GetPixel(u2, v0) +
    //          vq1 * GetPixel(u2, v1) +
    //          vq2 * GetPixel(u2, v2) +
    //          vq3 * GetPixel(u2, v3);
    //     D3 = vq0 * GetPixel(u3, v0) +
    //          vq1 * GetPixel(u3, v1) +
    //          vq2 * GetPixel(u3, v2) +
    //          vq3 * GetPixel(u3, v3);
    //     J.col(0) = dQ1(u) * D0 + dQ2(u) * D1 + dQ3(u) * D2 + dQ4(u) * D3;


    //     double uq0 = Q1(u);
    //     double uq1 = Q2(u);
    //     double uq2 = Q3(u);
    //     double uq3 = Q4(u);
    //     D0 = uq0 * GetPixel(u0, v0) +
    //          uq1 * GetPixel(u1, v0) +
    //          uq2 * GetPixel(u2, v0) +
    //          uq3 * GetPixel(u3, v0);
    //     D1 = uq0 * GetPixel(u0, v1) +
    //          uq1 * GetPixel(u1, v1) +
    //          uq2 * GetPixel(u2, v1) +
    //          uq3 * GetPixel(u3, v1);
    //     D2 = uq0 * GetPixel(u0, v2) +
    //          uq1 * GetPixel(u1, v2) +
    //          uq2 * GetPixel(u2, v2) +
    //          uq3 * GetPixel(u3, v2);
    //     D3 = uq0 * GetPixel(u0, v3) +
    //          uq1 * GetPixel(u1, v3) +
    //          uq2 * GetPixel(u2, v3) +
    //          uq3 * GetPixel(u3, v3);
    //     J.col(1) = dQ1(v) * D0 + dQ2(v) * D1 + dQ3(v) * D2 + dQ4(v) * D3;
    //     return J;
    // }
    // else
    //     throw std::runtime_error("Undefined interpolation scheme.");
}