/**
 * @file        segmentation.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#pragma once

#include <dfy/manifoldmesh.hpp>
#include <dfy/graph.hpp>
#include <dfy/region.hpp>


namespace dfy
{

typedef double (*MergeScore)(const dfy::ManifoldMesh&, const dfy::Region&, const dfy::Region&);
    
class Segmentation
{
private:
    const dfy::ManifoldMesh& m_Mesh;
    const dfy::Graph& m_Graph;
    std::vector<dfy::Region> m_Regions;
    Eigen::VectorXi m_Partitions;

public:
    Segmentation(const dfy::ManifoldMesh& M,
                 const dfy::Graph& G);
    Segmentation(const dfy::ManifoldMesh& M,
                 const dfy::Graph& G,
                 const Eigen::VectorXi& TriParts);
    ~Segmentation();

    void ComputeRegions(const Eigen::VectorXi& TriParts);
    const Eigen::VectorXi& GetTriParts() const;

    int NumRegions() const;
    const dfy::Region& GetRegion(int i) const;

    bool IsAllDisks() const;

    void MakeAllDisks(int SubSamples);
    void MergeRegions(dfy::MergeScore ScoreFun,
                      double Threshold);
    void CutToDisk(std::vector<int>& EdgeCuts,
                   const Eigen::VectorXd& EdgeWeights);
    void CutToDisk(std::vector<int>& EdgeCut);

    dfy::Graph DualGraph() const;
};

double MinPerimeterScore(const dfy::ManifoldMesh& M,
                         const dfy::Region& Ri,
                         const dfy::Region& Rj);

double MaxAvgDihedralAngle(const dfy::ManifoldMesh& M,
                           const dfy::Region& Ri,
                           const dfy::Region& Rj);

double MaxAvgCurvature(const dfy::ManifoldMesh& M,
                       const dfy::Region& Ri,
                       const dfy::Region& Rj);


} // namespace dfy
