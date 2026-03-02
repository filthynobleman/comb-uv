/**
 * @file        sampler.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#include <dfy/sampler.hpp>
#include <random>
#include <queue>


dfy::Sampler::Sampler(const dfy::Graph& G, int Seed) : m_G(G)
{
    std::mt19937 Eng(Seed);
    std::uniform_int_distribution<int> Dist(0, G.NumNodes() - 1);
    m_Samples.emplace_back(Dist(Eng));
    m_Partitions.setConstant(G.NumNodes(), 0);
    G.DijkstraDistance(m_Samples[0], m_Distances);
    m_Distances[m_Samples[0]] = -1;

    m_HDists = new dfy::MinHeap(m_Distances.data(), G.NumNodes(), true);
}

dfy::Sampler::Sampler(dfy::Sampler&& S) : m_G(S.m_G)
{
    m_Samples = std::move(S.m_Samples);
    m_Partitions = std::move(S.m_Partitions);
    m_Distances = std::move(S.m_Distances);
    m_HDists = S.m_HDists;
    S.m_HDists = nullptr;
}

dfy::Sampler::~Sampler()
{
    delete m_HDists;
}

double dfy::Sampler::GetDistance(int i) const { return m_Distances[i]; }
const Eigen::VectorXd &dfy::Sampler::GetDistances() const { return m_Distances; }
int dfy::Sampler::GetPartition(int i) const { return m_Partitions[i]; }
const Eigen::VectorXi &dfy::Sampler::GetPartitions() const { return m_Partitions; }

int dfy::Sampler::NumSamples() const { return m_Samples.size(); }
int dfy::Sampler::GetSample(int i) const { return m_Samples[i]; }
const std::vector<int> &dfy::Sampler::GetSamples() const { return m_Samples; }
int dfy::Sampler::FarthestVertex() const { return (int)(m_HDists->FindMin().second); }
void dfy::Sampler::AddSample() { AddSample(FarthestVertex()); }

void dfy::Sampler::AddSamples(int n)
{
    n += NumSamples();
    while (NumSamples() < n)
        AddSample();
}

void dfy::Sampler::AddSamples(const std::vector<int> &Samples)
{
    std::vector<int> S = Samples;
    std::sort(S.begin(), S.end());
    auto SEnd = std::unique(S.begin(), S.end());
    for (auto it = S.begin(); it != SEnd; it++)
        AddSample(*it);
}

void dfy::Sampler::AddSample(int NewSample)
{
    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> Q;
    // Let's be sure this specific sample will NEVER be kept again
    m_Distances[NewSample] = -1;
    m_Partitions[NewSample] = NumSamples();
    Q.emplace(m_Distances[NewSample], NewSample);
    while (!Q.empty())
    {
        std::pair<double, int> Next = Q.top();
        Q.pop();

        int Cur = Next.second;
        double W = Next.first;

        if (W < m_HDists->GetKey(Cur))
            m_HDists->SetKey(Cur, W);
        W = std::max(0.0, W);

        int Deg = m_G.NumAdjacents(Cur);
        for (int j = 0; j < Deg; ++j)
        {
            WEdge Neig = m_G.GetAdjacent(Cur, j);
            if (m_Distances[Neig.first] <= W + Neig.second)
                continue;
            m_Distances[Neig.first] = W + Neig.second;
            Q.emplace(m_Distances[Neig.first], Neig.first);
            m_Partitions[Neig.first] = NumSamples();
        }
    }

    m_Samples.emplace_back(NewSample);
}