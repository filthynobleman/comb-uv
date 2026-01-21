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

dfy::ARAPEmbedding::ARAPEmbedding(const dfy::Mesh& M) : dfy::TutteEmbedding(M), m_BFree(true) { }
dfy::ARAPEmbedding::ARAPEmbedding(const dfy::ARAPEmbedding& E) : dfy::TutteEmbedding(E), m_BFree(E.m_BFree) { }
dfy::ARAPEmbedding::ARAPEmbedding(dfy::ARAPEmbedding&& E) : dfy::TutteEmbedding(E), m_BFree(E.m_BFree) { }
dfy::ARAPEmbedding::~ARAPEmbedding() { }

void dfy::ARAPEmbedding::MapBoundary(dfy::BoundaryMap BMap)
{
    m_BFree = false;
    dfy::TutteEmbedding::MapBoundary(BMap);
}

bool dfy::ARAPEmbedding::Compute()
{
    if (m_BFree)
        dfy::TutteEmbedding::MapBoundary();
    if (!dfy::TutteEmbedding::Compute())
        return false;

    igl::ARAPData data;
    data.with_dynamics = !m_BFree;
    data.max_iter = 100;
    Eigen::VectorXi b(0);
    Eigen::MatrixXd bc(0, 0);
    if (!m_BFree)
    {
        b = Eigen::VectorXi::Map(GetBoundaryLoop().data(), GetBoundaryLoop().size());
        bc = GetUV()(b, Eigen::all);
        data.max_iter = 1000;
        data.ym = 0.1;
    }
    igl::arap_precomputation(GetMesh().Vertices(), GetMesh().Triangles(), 2, b, data);
    igl::arap_solve(bc, data, GetUV());

    if (m_BFree)
    {
        GetUV().col(0) = GetUV().col(0).array() - GetUV().col(0).minCoeff();
        GetUV().col(1) = GetUV().col(1).array() - GetUV().col(1).minCoeff();
        GetUV() = GetUV().array() / GetUV().maxCoeff();
    }
    return true;
}