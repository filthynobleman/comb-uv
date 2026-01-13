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
public:
    ARAPEmbedding(const dfy::Mesh& M);
    ARAPEmbedding(const dfy::ARAPEmbedding& E);
    ARAPEmbedding(dfy::ARAPEmbedding&& E);
    ~ARAPEmbedding();

    virtual bool Compute() override;
};


} // namespace dfy
