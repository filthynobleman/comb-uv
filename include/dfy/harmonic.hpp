/**
 * @file        harmonic.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#pragma once

#include <dfy/embedding.hpp>



namespace dfy
{
    
class HarmonicEmbedding : public dfy::Embedding
{
public:
    HarmonicEmbedding(const dfy::Mesh& M);
    HarmonicEmbedding(const dfy::HarmonicEmbedding& E);
    HarmonicEmbedding(dfy::HarmonicEmbedding&& E);
    ~HarmonicEmbedding();

    virtual bool Compute() override;
};


} // namespace dfy
