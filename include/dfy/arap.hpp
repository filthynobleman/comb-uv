/**
 * @file        arap.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#pragma once

#include <dfy/tutte.hpp>



namespace dfy
{
    
class ARAPEmbedding : public dfy::TutteEmbedding
{
private:
    bool m_BFree;
public:
    ARAPEmbedding(const dfy::Mesh& M);
    ARAPEmbedding(const dfy::ARAPEmbedding& E);
    ARAPEmbedding(dfy::ARAPEmbedding&& E);
    ~ARAPEmbedding();

    void SetBoundaryFree();

    virtual void MapBoundary(dfy::BoundaryMap BMap = dfy::BoundaryMap::CIRCLE) override;

    virtual bool Compute() override;
};


} // namespace dfy
