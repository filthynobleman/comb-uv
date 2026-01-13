/**
 * @file        tutte.hpp
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
    
class TutteEmbedding : public dfy::Embedding
{
public:
    TutteEmbedding(const dfy::Mesh& M);
    TutteEmbedding(const dfy::TutteEmbedding& E);
    TutteEmbedding(dfy::TutteEmbedding&& E);
    ~TutteEmbedding();

    virtual bool Compute() override;
};


} // namespace dfy
