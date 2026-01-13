/**
 * @file        conformal.hpp
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
    
class ConformalEmbedding : public dfy::Embedding
{
public:
    ConformalEmbedding(const dfy::Mesh& M);
    ConformalEmbedding(const dfy::ConformalEmbedding& E);
    ConformalEmbedding(dfy::ConformalEmbedding&& E);
    ~ConformalEmbedding();

    virtual bool Compute() override;
};


} // namespace dfy
