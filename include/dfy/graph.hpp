/**
 * @file        graph.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#pragma once

#include <dfy/utils.hpp>
#include <dfy/mesh.hpp>
#include <dfy/manifoldmesh.hpp>

namespace dfy
{
    
/**
 * @brief       Weighted edge of a graph.
 * 
 * @details     A weighted edge of a graph is defined as a destintation (int) and a cost (double).
 */
typedef std::pair<int, double> WEdge;

/**
 * @brief       Path inside a weighted graph.
 * 
 * @details     A path inside a weighted graph is defined as sequence of weighted edges
 *              and the sum of the total weights. The starting node has cost zero.
 */
typedef std::pair<double, std::vector<dfy::WEdge>> Path;


class Graph
{
private:
    std::vector<int> m_Idxs;
    std::vector<dfy::WEdge> m_Adjs;
    Graph();

public:
    Graph(const std::vector<std::pair<int, int>>& Edges,
          const std::vector<double>& Weights = {});
    Graph(const dfy::Graph& G);
    Graph(dfy::Graph&& G);
    dfy::Graph& operator=(const dfy::Graph& G);
    dfy::Graph& operator=(dfy::Graph&& G);
    ~Graph();


    int NumNodes() const;
    int NumAdjacents(int i) const;
    int NumEdges() const;

    void CollapseEdge(int i, int j);

    const dfy::WEdge& GetAdjacent(int node_i, int adj_i) const;
    dfy::Path DijkstraPath(int src, int dst) const;
    Eigen::VectorXd DijkstraDistance(int src) const;
    void DijkstraDistance(int src, Eigen::VectorXd& Distances) const;
    Eigen::VectorXd DijkstraDistance(int src, const Eigen::VectorXi& Tag, int Filter) const;
    std::vector<int> ConnectedComponents() const;
    std::vector<std::pair<int, int>> MaxSpanTree() const;

    dfy::Graph SubGraph(const std::vector<int>& Indices) const;
};



typedef double (*M2GDist)(const dfy::Mesh&, int, int);
typedef double (*DM2GDist)(const dfy::ManifoldMesh&, int, int);

double EuclideanDistance(const dfy::Mesh& M, int i, int j);
double DualEuclideanDistance(const dfy::ManifoldMesh& M, int i, int j);
double GeodesicDistance(const dfy::ManifoldMesh& M, int i, int j);
double AngularDistance(const dfy::Mesh& M, int i, int j);
double DualAngularDistance(const dfy::ManifoldMesh& M, int i, int j);
double DualCurvatureDistance(const dfy::ManifoldMesh& M, int i, int j);
double ConstantDistance(const dfy::Mesh&, int i, int j);
double DualConstantDistance(const dfy::ManifoldMesh&, int i, int j);


dfy::Graph MeshToGraph(const dfy::Mesh& M,
                       M2GDist Dist = EuclideanDistance);
dfy::Graph DualMeshToGraph(const dfy::ManifoldMesh& M,
                           DM2GDist Dist = DualEuclideanDistance);


} // namespace dfy
