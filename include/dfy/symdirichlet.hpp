/**
 * @file        symdirichlet.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#pragma once

#include <dfy/harmonic.hpp>



namespace dfy
{
    
class SymDirichletEmbedding : public dfy::HarmonicEmbedding
{
private:
    bool m_BFree;
public:
    SymDirichletEmbedding(const dfy::Mesh& M);
    SymDirichletEmbedding(const dfy::SymDirichletEmbedding& E);
    SymDirichletEmbedding(dfy::SymDirichletEmbedding&& E);
    ~SymDirichletEmbedding();

    void SetBoundaryFree();

    virtual void MapBoundary(dfy::BoundaryMap BMap = dfy::BoundaryMap::CIRCLE) override;

    virtual bool Compute() override;
};


} // namespace dfy
