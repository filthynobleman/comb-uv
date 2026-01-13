/**
 * @file        harmonic.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#include <dfy/harmonic.hpp>

#include <igl/harmonic.h>

dfy::HarmonicEmbedding::HarmonicEmbedding(const dfy::Mesh& M) : dfy::Embedding(M) { }
dfy::HarmonicEmbedding::HarmonicEmbedding(const dfy::HarmonicEmbedding& E) : dfy::Embedding(E) { }
dfy::HarmonicEmbedding::HarmonicEmbedding(dfy::HarmonicEmbedding&& E) : dfy::Embedding(E) { }
dfy::HarmonicEmbedding::~HarmonicEmbedding() { }

bool dfy::HarmonicEmbedding::Compute()
{
    Eigen::VectorXi BLoop = Eigen::VectorXi::Map(GetBoundaryLoop().data(), GetBoundaryLoop().size());
    if (!igl::harmonic(GetMesh().Vertices(), GetMesh().Triangles(), BLoop, GetUV()(BLoop, Eigen::all).eval(), 1, GetUV()))
        return false;
    for (int i = 0; i < GetUV().rows(); ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            if (std::isinf(GetUV()(i, j)) || std::isnan(GetUV()(i, j)))
                return false;
        }
    }
    return true;
}