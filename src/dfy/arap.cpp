/**
 * @file        arap.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#include <dfy/arap.hpp>

#include <igl/arap.h>

dfy::ARAPEmbedding::ARAPEmbedding(const dfy::Mesh& M) : dfy::TutteEmbedding(M) { }
dfy::ARAPEmbedding::ARAPEmbedding(const dfy::ARAPEmbedding& E) : dfy::TutteEmbedding(E) { }
dfy::ARAPEmbedding::ARAPEmbedding(dfy::ARAPEmbedding&& E) : dfy::TutteEmbedding(E) { }
dfy::ARAPEmbedding::~ARAPEmbedding() { }

bool dfy::ARAPEmbedding::Compute()
{
    if (!dfy::TutteEmbedding::Compute())
        return false;

    igl::ARAPData data;
    data.with_dynamics = true;
    data.max_iter = 100;
    Eigen::VectorXi b(0);
    Eigen::MatrixXd bc(0, 0);
    igl::arap_precomputation(GetMesh().Vertices(), GetMesh().Triangles(), 2, b, data);
    igl::arap_solve(bc, data, GetUV());

    GetUV().col(0) = GetUV().col(0).array() - GetUV().col(0).minCoeff();
    GetUV().col(1) = GetUV().col(1).array() - GetUV().col(1).minCoeff();
    GetUV() = GetUV().array() / GetUV().maxCoeff();
    return true;
}