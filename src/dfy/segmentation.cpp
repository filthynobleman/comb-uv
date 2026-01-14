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

dfy::Graph dfy::Segmentation::DualGraph() const
{
    std::vector<std::pair<int, int>> Edges;
    std::vector<double> Weights;
    std::vector<int> Eij;
    for (int i = 0; i < NumRegions(); ++i)
    {
        const auto& Ri = GetRegion(i);
        for (int j = i + 1; j < NumRegions(); ++j)
        {
            const auto& Rj = GetRegion(j);
            auto Rij = dfy::LineIntersection(Ri, Rj);
            if (Rij.IsEmpty())
                continue;
            Rij.Edges(Eij);
            double w = m_Mesh.EdgeLengths()(Eij).sum();
            Edges.emplace_back(i, j);
            Weights.emplace_back(w);
        }
    }

    return dfy::Graph(Edges, Weights);
}

void dfy::Segmentation::CutToDisk(std::vector<int> &EdgeCut)
{
    EdgeCut.clear();

    std::vector<std::pair<int, int>> RegEdges;
    RegEdges.reserve(6 * NumRegions());
    std::vector<double> Weights;
    Weights.reserve(6 * NumRegions());

    // Compute the length of the mesh edges
    const Eigen::VectorXd& ELens = m_Mesh.EdgeLengths();
    // std::cout << ELens.rows() << std::endl;

    // Compute the adjacency graph
    std::vector<int> Eij;
    std::vector<int> Vij;
    Eigen::VectorXi ECut;
    ECut.setZero(m_Mesh.NumEdges());
    for (int i = 0; i < NumRegions(); ++i)
    {
        for (int j = i + 1; j < NumRegions(); ++j)
        {
            auto Rij = dfy::LineIntersection(GetRegion(i), GetRegion(j));
            // Empty intersection
            if (Rij.EulerCharacteristic() == 0)
                continue;
            // Add to list of potential cuts
            Rij.Edges(Eij);
            ECut(Eij).setOnes();
            // std::cout << Rij.EulerCharacteristic() << std::endl;
            // Single segment is easy
            if (Rij.EulerCharacteristic() == 1)
            {
                RegEdges.emplace_back(i, j);
                // Negated weight, as we want the maximum spanning tree
                Weights.emplace_back(ELens(Eij).sum());
                continue;
            }
            // With multiple segments, we first need the connected components
            // Build vertex-edge adjacency; each vertex has AT MOST 2 neighboring edges
            Rij.Vertices(Vij);
            std::map<int, std::pair<int, int>> V2E;
            std::pair<int, int> MinusOne{ -1, -1 };
            for (int k = 0; k < Eij.size(); ++k)
            {
                int v = m_Mesh.Edges()(Eij[k], 0);
                if (V2E.find(v) == V2E.end())
                    V2E.emplace(v, MinusOne);
                if (V2E[v].first == -1)
                    V2E[v].first = k;
                else
                    V2E[v].second = k;
                
                v = m_Mesh.Edges()(Eij[k], 1);
                if (V2E.find(v) == V2E.end())
                    V2E.emplace(v, MinusOne);
                if (V2E[v].first == -1)
                    V2E[v].first = k;
                else
                    V2E[v].second = k;
            }
            // std::cout << "V2E built" << std::endl;
            // Build edge-edge graph; each edge has AT MOST 2 neighboring edges
            std::vector<std::pair<int, int>> E2E;
            for (auto it : V2E)
            {
                if (it.second.first == -1 || it.second.second == -1)
                    continue;
                E2E.push_back(it.second);
            }
            // If all the boundary components are one-edge long, we just cut through all of them
            if (E2E.empty())
            {
                EdgeCut.insert(EdgeCut.end(), Eij.begin(), Eij.end());
                continue;
            }
            // std::cout << "E2E built" << std::endl;
            dfy::Graph EEG(E2E);
            // std::cout << "Graph built" << std::endl;
            // Find connected components and partition the edges
            std::vector<int> ECC = EEG.ConnectedComponents();
            auto it = std::max_element(ECC.begin(), ECC.end());
            // std::cout << ECC.size() << " == " << Eij.size() << std::endl;
            // std::cout << (*it) << " == " << Rij.EulerCharacteristic() << std::endl;
            std::vector<std::vector<int>> Chains;
            Chains.resize(std::max(Rij.EulerCharacteristic(), (*it) + 1));
            for (int k = 0; k < ECC.size(); ++k)
            {
                
                Chains[ECC[k]].emplace_back(Eij[k]);
            }
            // std::cout << "Chains built" << std::endl;
            // Find longest boundary
            int MaxChainIdx = 0;
            double MaxChainLen = ELens(Chains[0]).sum();
            for (int k = 1; k < Chains.size(); ++k)
            {
                double ChainLen = ELens(Chains[k]).sum();
                if (ChainLen <= MaxChainLen)
                    continue;
                MaxChainLen = ChainLen;
                MaxChainIdx = k;
            }
            // std::cout << "Max chain found" << std::endl;
            // Add longest boundary to graph
            RegEdges.emplace_back(i, j);
            Weights.emplace_back(MaxChainLen);
            // All other edges are to cut for sure
            for (int k = 0; k < Chains.size(); ++k)
            {
                if (k == MaxChainIdx)
                    continue;
                EdgeCut.insert(EdgeCut.end(), Chains[k].begin(), Chains[k].end());
            }
            // std::cout << "Added edges to cut" << std::endl;
        }
    }
    dfy::Graph AdjGraph(RegEdges, Weights);

    // Get edges in the minimum spanning tree
    std::vector<std::pair<int, int>> MST = AdjGraph.MinSpanTree();

    // Remove MST edges from potential cuts
    for (auto e : MST)
    {
        auto Rij = dfy::LineIntersection(GetRegion(e.first), GetRegion(e.second));
        std::vector<int> Eij;
        Rij.Edges(Eij);
        ECut(Eij).setZero();
    }

    // Add potential cuts to cut list
    for (int e = 0; e < m_Mesh.NumEdges(); ++e)
    {
        if (ECut[e] > 0)
            EdgeCut.emplace_back(e);
    }
}