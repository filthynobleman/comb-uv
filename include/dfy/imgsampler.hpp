/**
 * @file        imgsampler.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-15
 */
#pragma once

#include <dfy/utils.hpp>
#include <dfy/gimage.hpp>


namespace dfy
{

enum ImgInterpType
{
    LINEAR,
    CUBIC
};
    
class ImageSampler
{
private:
    int m_Width;
    int m_Height;
    Eigen::Vector3d* m_Data;

public:
    ImageSampler(const std::string& Filename);
    ImageSampler(const dfy::GImage& Img);
    ImageSampler(const dfy::ImageSampler& IS);
    ImageSampler(dfy::ImageSampler&& IS);
    dfy::ImageSampler& operator=(const dfy::ImageSampler& IS);
    dfy::ImageSampler& operator=(dfy::ImageSampler&& IS);
    ~ImageSampler();

    int GetWidth() const;
    int GetHeight() const;

    const Eigen::Vector3d& GetPixel(int i, int j) const;
    Eigen::Vector3d Sample(double u, double v,
                           dfy::ImgInterpType InterpType = dfy::LINEAR) const;
    Eigen::Vector3d Sample(const Eigen::Vector2d& uv,
                           dfy::ImgInterpType InterpType = dfy::LINEAR) const;
    Eigen::Matrix<double, 3, 2> Jacobian(double u, double v,
                                         dfy::ImgInterpType InterpType = dfy::LINEAR) const;
    Eigen::Matrix<double, 3, 2> Jacobian(const Eigen::Vector2d& uv,
                                         dfy::ImgInterpType InterpType = dfy::LINEAR) const;
};


} // namespace dfy
