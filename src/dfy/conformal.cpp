/**
 * @file        conformal.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#include <dfy/conformal.hpp>

#include <igl/lscm.h>

dfy::ConformalEmbedding::ConformalEmbedding(const dfy::Mesh& M) : dfy::Embedding(M) { }
dfy::ConformalEmbedding::ConformalEmbedding(const dfy::ConformalEmbedding& E) : dfy::Embedding(E) { }
dfy::ConformalEmbedding::ConformalEmbedding(dfy::ConformalEmbedding&& E) : dfy::Embedding(E) { }
dfy::ConformalEmbedding::~ConformalEmbedding() { }

bool dfy::ConformalEmbedding::Compute()
{
    if (!igl::lscm(GetMesh().Vertices(), GetMesh().Triangles(), GetUV()))
        return false;

    GetUV().col(0) = GetUV().col(0).array() - GetUV().col(0).minCoeff();
    GetUV().col(1) = GetUV().col(1).array() - GetUV().col(1).minCoeff();
    GetUV() = GetUV().array() / GetUV().maxCoeff();

    return true;
}