/**
 * @file        embedding.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#pragma once


#include <dfy/utils.hpp>
#include <dfy/mesh.hpp>


namespace dfy
{

enum BoundaryMap
{
    CIRCLE,
    SQUARE,
    TRIANGLE,
    HEXAGON
};
    
class Embedding
{
private:
    const dfy::Mesh& m_Mesh;

    std::vector<int> m_BLoop;
    Eigen::MatrixXd m_UV;

protected:
    const dfy::Mesh& GetMesh() const;
    std::vector<int>& GetBoundaryLoop();
    Eigen::MatrixXd& GetUV();


public:
    Embedding(const dfy::Mesh& M);
    Embedding(const dfy::Embedding& E);
    Embedding(dfy::Embedding&& E);
    ~Embedding();
    
    const std::vector<int>& GetBoundaryLoop() const;
    const Eigen::MatrixXd& GetUV() const;

    void MapBoundary(dfy::BoundaryMap BMap = dfy::BoundaryMap::CIRCLE);
    virtual bool Compute() = 0;
};

} // namespace dfy
