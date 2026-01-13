/**
 * @file        segmentation.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#include <dfy/segmentation.hpp>

#include <dfy/utils.hpp>
#include <dfy/sampler.hpp>


dfy::Segmentation::Segmentation(const dfy::ManifoldMesh &M,
                                const dfy::Graph& G)
    : m_Mesh(M), m_Graph(G)
{
    ComputeRegions(Eigen::VectorXi::Zero(m_Mesh.NumTriangles()).eval());
}

dfy::Segmentation::Segmentation(const dfy::ManifoldMesh &M,
                                const dfy::Graph& G,
                                const Eigen::VectorXi &TriParts)
    : m_Mesh(M), m_Graph(G)
{
    ComputeRegions(TriParts);
}

dfy::Segmentation::~Segmentation() { }

void dfy::Segmentation::ComputeRegions(const Eigen::VectorXi &TriParts)
{
    m_Partitions = TriParts;
    m_Regions.clear();
    m_Regions.resize(m_Partitions.maxCoeff() + 1);

    for (int i = 0; i < m_Mesh.NumTriangles(); ++i)
    {
        m_Regions[m_Partitions[i]].AddFace(i);
        
        for (int j = 0; j < 3; ++j)
        {
            if (m_Mesh.TriEdgeAdj()(i, j) >= 0)
                m_Regions[m_Partitions[i]].AddEdge(m_Mesh.TriEdgeAdj()(i, j));
            m_Regions[m_Partitions[i]].AddVertex(m_Mesh.Triangles()(i, j));
        }
    }
}

const Eigen::VectorXi &dfy::Segmentation::GetTriParts() const { return m_Partitions; }
int dfy::Segmentation::NumRegions() const { return m_Regions.size(); }
const dfy::Region &dfy::Segmentation::GetRegion(int i) const { return m_Regions[i]; }

bool dfy::Segmentation::IsAllDisks() const
{
    for (const auto& r : m_Regions)
    {
        if (!r.IsDisk())
            return false;
    }
    return true;
}

void dfy::Segmentation::MakeAllDisks(int SubSamples)
{
    bool Changes = true;
    while (Changes)
    {
        Changes = false;
        Eigen::VectorXi NewParts = m_Partitions;
        int NewNumRegs = NumRegions();
        for (int rid = 0; rid < NumRegions(); ++rid)
        {
            const auto& r = GetRegion(rid);
            if (r.IsDisk())
                continue;
            Changes = true;
            
            std::vector<int> RFaces;
            r.Faces(RFaces);
            dfy::Graph SG = m_Graph.SubGraph(RFaces);
            dfy::Sampler Smpl(SG);
            while (Smpl.NumSamples() < std::min(SubSamples, SG.NumNodes()))
                Smpl.AddSample();
            for (int i = 0; i < RFaces.size(); ++i)
            {
                if (Smpl.GetPartition(i) == 0)
                    continue;
                NewParts[RFaces[i]] = NewNumRegs + Smpl.GetPartition(i) - 1;
            }
            NewNumRegs += Smpl.NumSamples() - 1;
        }

        if (Changes)
            ComputeRegions(NewParts);
    }
}