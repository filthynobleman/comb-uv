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

#include <queue>
#include <stack>

#include <iostream>
#include <dfy/io.hpp>


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
        // std::cout << NumRegions() << std::endl;
        Changes = false;
        Eigen::VectorXi NewParts = m_Partitions;
        int NewNumRegs = NumRegions();
        int OldRegions = NumRegions();
        // int ExpNewRegs = 0;
        // for (int rid = 0; rid < NumRegions(); ++rid)
        //     std::cout << rid << ": " << GetRegion(rid).NumFaces() << std::endl;
        // std::cout << std::endl;
        // m_Regions.reserve(NewNumRegs + ExpNewRegs);
        for (int rid = 0; rid < NumRegions(); ++rid)
        {
            const auto& r = GetRegion(rid);
            if (r.IsDisk())
                continue;
            Changes = true;
            
            std::vector<int> RFaces;
            r.Faces(RFaces);
            dfy::Graph SG = m_Graph.SubGraph(RFaces);
            // First, check connected components are not creating issues
            std::vector<int> SGCC = SG.ConnectedComponents();
            int NCC = 0;
            for (int cc : SGCC)
                NCC = std::max(cc + 1, NCC);
            // If we have multiple connected components, just separate
            int OldNumRegs = NewNumRegs;
            if (NCC > 1)
            {
                for (int i = 0; i < RFaces.size(); ++i)
                {
                    if (SGCC[i] == 0)
                        NewParts[RFaces[i]] = rid;
                    else
                        NewParts[RFaces[i]] = NewNumRegs + SGCC[i] - 1;
                }
                NewNumRegs += NCC - 1;
            }
            // Otherwise, break by sampling
            else
            {
                dfy::Sampler Smpl(SG);
                while (Smpl.NumSamples() < std::min(SubSamples, SG.NumNodes()))
                    Smpl.AddSample();
                for (int i = 0; i < RFaces.size(); ++i)
                {
                    if (Smpl.GetPartition(i) == 0)
                        NewParts[RFaces[i]] = rid;
                    else
                        NewParts[RFaces[i]] = NewNumRegs + Smpl.GetPartition(i) - 1;
                }
                NewNumRegs += Smpl.NumSamples() - 1;
            }
            // Update regions data structure
            m_Regions[rid].Clear();
            // Preallocation and tricky memory management saves about 30-40% of computation
            int rV = 2 * r.NumVertices() / SubSamples;
            int rE = 2 * r.NumEdges() / SubSamples;
            int rT = 2 * r.NumFaces() / SubSamples;
            m_Regions[rid].Reserve(rV, rE, rT);
            if (m_Regions.capacity() < NewNumRegs)
                m_Regions.reserve(NewNumRegs + 10 * SubSamples);
            m_Regions.resize(NewNumRegs);
            for (int i = OldNumRegs; i < NewNumRegs; ++i)
                m_Regions[i].Reserve(rV, rE, rT);
            for (int i = 0; i < RFaces.size(); ++i)
            {
                int Part = NewParts[RFaces[i]];
                m_Regions[Part].AddFace(RFaces[i]);
                for (int j = 0; j < 3; ++j)
                {
                    // if (m_Mesh.TriEdgeAdj()(RFaces[i], j) >= 0)
                    m_Regions[Part].AddEdge(m_Mesh.TriEdgeAdj()(RFaces[i], j));
                    m_Regions[Part].AddVertex(m_Mesh.Triangles()(RFaces[i], j));
                }
            }
        }

        m_Partitions = NewParts;


        // if (Changes)
        //     ComputeRegions(NewParts);
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

void dfy::Segmentation::MergeRegions(dfy::MergeScore ScoreFun, 
                                     double Threshold)
{
    dfy::Graph G = DualGraph();

    bool SomethingMerged = false;
    do
    {
        // Scored approach:
        // Assign a score to each edge, and select the best scoring edge
        // Score must not be less than a threshold
        // Score = Boundary length / maximum area
        SomethingMerged = false;
        std::priority_queue<std::pair<double, std::pair<int, int>>> Q;
        for (int i = 0; i < G.NumNodes(); ++i)
        {
            int nadj = G.NumAdjacents(i);
            for (int jj = 0; jj < nadj; ++jj)
            {
                int j = G.GetAdjacent(i, jj).first;
                if (j < i)
                    continue;
                dfy::Region Intersection = dfy::LineIntersection(m_Regions[i], m_Regions[j]);
                if (Intersection.NumEdges() == 0)
                    continue;
                if (Intersection.EulerCharacteristic() != 1)
                    continue;
                double Score = ScoreFun(m_Mesh, m_Regions[i], m_Regions[j]);
                if (Score < Threshold)
                {
                    // std::cout << Score << std::endl;
                    continue;
                }
                std::pair<int, int> ij{ i, j };
                Q.emplace(Score, ij);
            }
        }
        if (Q.empty())
            break;
        
        SomethingMerged = true;
        double Score;
        std::pair<int, int> e;
        std::tie(Score, e) = Q.top();
        int i, j;
        std::tie(i, j) = e;
        dfy::Region Union = dfy::Union(m_Regions[i], m_Regions[j]);
        if (Union.EulerCharacteristic() != 1)
            continue;

        m_Regions[i] = Union;
        m_Regions.erase(m_Regions.begin() + j);

        for (int f = 0; f < m_Partitions.rows(); ++f)
        {
            if (m_Partitions[f] == j)
                m_Partitions[f] = i;
            else if (m_Partitions[f] >= j)
                m_Partitions[f] -= 1;
        }


        SomethingMerged = true;
        try
        {
            G.CollapseEdge(i, j);
        }
        catch (std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    } while(SomethingMerged);
}

void dfy::Segmentation::CutToDisk(std::vector<int> &EdgeCut)
{
    CutToDisk(EdgeCut, m_Mesh.EdgeLengths());
}

void dfy::Segmentation::CutToDisk(std::vector<int> &EdgeCut,
                                  const Eigen::VectorXd& EdgeWeights)
{
    EdgeCut.clear();

    std::vector<std::pair<int, int>> RegEdges;
    RegEdges.reserve(6 * NumRegions());
    std::vector<double> Weights;
    Weights.reserve(6 * NumRegions());

    // Compute the length of the mesh edges
    const Eigen::VectorXd& ELens = EdgeWeights;
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
            ECut(Eij).setConstant(2);
            // Single segment is easy
            if (Rij.EulerCharacteristic() == 1)
            {
                RegEdges.emplace_back(i, j);
                Weights.emplace_back(ELens(Eij).sum());
                ECut(Eij).setOnes();
                continue;
            }
            // With multiple segments, we first need the connected components
            // Build vertex-edge adjacency
            // We are on a boundary component (a line), so each vertex has AT MOST 2 neighboring edges
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
            // Build edge-edge graph
            // We are on a boundary component (a line), so each edge has AT MOST 2 neighboring edges
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
                Chains[ECC[k]].emplace_back(Eij[k]);
            // Find longest boundary
            int MaxChainIdx = 0;
            double MaxChainLen = ELens(Chains[0]).sum();
            for (int k = 0; k < Chains.size(); ++k)
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
            ECut(Chains[MaxChainIdx]).setOnes();
            // // All other edges are to cut for sure
            // for (int k = 0; k < Chains.size(); ++k)
            // {
            //     if (k == MaxChainIdx)
            //         continue;
            //     EdgeCut.insert(EdgeCut.end(), Chains[k].begin(), Chains[k].end());
            // }
            // std::cout << "Added edges to cut" << std::endl;
        }
    }
    dfy::Graph AdjGraph(RegEdges, Weights);

    // Get edges in the minimum spanning tree
    std::vector<std::pair<int, int>> MST = AdjGraph.MaxSpanTree();

    // Remove MST edges from potential cuts
    for (auto e : MST)
    {
        auto Rij = dfy::LineIntersection(GetRegion(e.first), GetRegion(e.second));
        std::vector<int> Eij;
        Rij.Edges(Eij);
        ECut(Eij) = ECut(Eij).array() - 1;
    }

    // Add potential cuts to cut list
    for (int e = 0; e < m_Mesh.NumEdges(); ++e)
    {
        if (ECut[e] > 0)
            EdgeCut.emplace_back(e);
    }
}

double dfy::MinPerimeterScore(const dfy::ManifoldMesh &M, 
                              const dfy::Region &Ri, 
                              const dfy::Region &Rj)
{
    std::vector<int> FID;
    Ri.Faces(FID);
    double Ai = M.FaceAreas()(FID).sum();
    Rj.Faces(FID);
    double Aj = M.FaceAreas()(FID).sum();
    auto Rij = dfy::LineIntersection(Ri, Rj);
    Rij.Edges(FID);
    double Eij = M.EdgeLengths()(FID).sum();

    return Eij / std::sqrt(std::max(Ai, Aj));
}


#define _USE_MATH_DEFINES
#include <math.h>

double dfy::MaxAvgDihedralAngle(const dfy::ManifoldMesh &M, 
                                const dfy::Region &Ri, 
                                const dfy::Region &Rj)
{
    std::vector<int> EID;
    auto Rij = dfy::LineIntersection(Ri, Rj);
    Rij.Edges(EID);
    double TotDihedral = 0.0;
    for (int e : EID)
    {
        double Dot = M.FaceNormals().row(M.EdgeTriAdj()(e, 0)).dot(M.FaceNormals().row(M.EdgeTriAdj()(e, 1)));
        // Clamp to [-1, 1]
        Dot = std::min(1.0 - 1e-12, std::max(-1.0 + 1e-12, Dot));
        TotDihedral += std::acos(Dot);
    }
    if (TotDihedral < 1e-9)
        return 0.0;
    if (EID.size() < 1e-9)
        return 0.0;
    // std::cout << (M_PI - TotDihedral / EID.size()) / M_PI * 180 << std::endl;
    return (M_PI - TotDihedral / EID.size()) / M_PI * 180.0;
}

double dfy::MaxAvgCurvature(const dfy::ManifoldMesh &M, 
                            const dfy::Region &Ri, 
                            const dfy::Region &Rj)
{
    std::vector<int> EID;
    auto Rij = dfy::LineIntersection(Ri, Rj);
    Rij.Edges(EID);
    std::vector<int> VID;
    VID.reserve(2 * EID.size());
    for (int e : EID)
    {
        VID.emplace_back(M.Edges()(e, 0));
        VID.emplace_back(M.Edges()(e, 1));
    }
    std::sort(VID.begin(), VID.end());
    auto VEnd = std::unique(VID.begin(), VID.end());
    VID.erase(VEnd, VID.end());
    double GC = (M.GaussianCurvature()(VID).cwiseAbs() + M.MeanCurvature()(VID).cwiseAbs()).mean();
    return 1.0 / GC;
}