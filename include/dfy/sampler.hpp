/**
 * @file        voronoifps.hpp
 * 
 * @brief       Declaration of class dfy::Sampler.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2024-01-15
 */
#pragma once

#include <dfy/graph.hpp>
#include <dfy/mesh.hpp>
#include <dfy/minheap.hpp>


namespace dfy
{
    
class Sampler
{
private:
    const dfy::Graph& m_G;
    std::vector<int> m_Samples;
    Eigen::VectorXi m_Partitions;
    Eigen::VectorXd m_Distances;
    dfy::MinHeap* m_HDists;


public:
    Sampler(const dfy::Graph& G,
            int Seed = 0);
    Sampler(dfy::Sampler&& VP);
    ~Sampler();

    double GetDistance(int i) const;
    const Eigen::VectorXd& GetDistances() const;
    int GetPartition(int i) const;
    const Eigen::VectorXi& GetPartitions() const;
    
    int NumSamples() const;
    int GetSample(int i) const;
    const std::vector<int>& GetSamples() const;

    int FarthestVertex() const;
    void AddSample();
    void AddSample(int NewSample);
    void AddSamples(int n);
    void AddSamples(const std::vector<int>& Samples);
};



} // namespace dfy