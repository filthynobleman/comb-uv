/**
 * @file        symdirichlet.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#include <dfy/symdirichlet.hpp>
#include <igl/slim.h>

dfy::SymDirichletEmbedding::SymDirichletEmbedding(const dfy::Mesh& M) : dfy::HarmonicEmbedding(M), m_BFree(true) { }
dfy::SymDirichletEmbedding::SymDirichletEmbedding(const dfy::SymDirichletEmbedding& E) : dfy::HarmonicEmbedding(E), m_BFree(E.m_BFree) { }
dfy::SymDirichletEmbedding::SymDirichletEmbedding(dfy::SymDirichletEmbedding&& E) : dfy::HarmonicEmbedding(E), m_BFree(E.m_BFree) { }
dfy::SymDirichletEmbedding::~SymDirichletEmbedding() { }

void dfy::SymDirichletEmbedding::MapBoundary(dfy::BoundaryMap BMap)
{
    m_BFree = false;
    dfy::HarmonicEmbedding::MapBoundary(BMap);
}

bool dfy::SymDirichletEmbedding::Compute()
{
    if (m_BFree)
        dfy::HarmonicEmbedding::MapBoundary();
    if (!dfy::HarmonicEmbedding::Compute())
        return false;

    Eigen::VectorXi BL = Eigen::VectorXi::Map(GetBoundaryLoop().data(), GetBoundaryLoop().size());
    Eigen::MatrixXd BC = GetUV()(GetBoundaryLoop(), Eigen::all);
    igl::SLIMData SData;
    GetUV().setConstant(GetMesh().NumVertices(), 2, 0.5);
    GetUV()(BL, Eigen::all) = BC;
    igl::slim_precompute(GetMesh().Vertices(), GetMesh().Triangles(),
                         GetUV(), SData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET,
                         BL, BC, 0);

    igl::slim_solve(SData, 10);
    Eigen::MatrixXd Tmp = SData.V_o;
    std::cout << SData.energy << std::endl;
    std::cout << (Tmp - GetUV()).eval().norm() << std::endl;
    GetUV() = Tmp;

    if (m_BFree)
    {
        GetUV().col(0) = GetUV().col(0).array() - GetUV().col(0).minCoeff();
        GetUV().col(1) = GetUV().col(1).array() - GetUV().col(1).minCoeff();
        GetUV() = GetUV().array() / GetUV().maxCoeff();
    }
    return true;
}