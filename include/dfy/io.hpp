/**
 * @file        io.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#pragma once

#include <dfy/utils.hpp>
#include <dfy/graph.hpp>
#include <dfy/mesh.hpp>

#include <ostream>
#include <fstream>


namespace dfy
{
    
bool ExportMesh(const std::string& Filename,
                const dfy::Mesh& M);

bool ExportMesh(const std::string& Filename,
                const dfy::Mesh& M,
                const Eigen::MatrixXd& UV,
                const Eigen::MatrixXi& TUV,
                bool Smooth = true);

bool ExportPointCloud(const std::string& Filename,
                      const Eigen::MatrixXd& Points);

bool ExportGraph(const std::string& Filename,
                 const dfy::Graph& G,
                 const Eigen::MatrixXd& V);

bool ExportCut(const std::string& Filename,
               const dfy::ManifoldMesh& M,
               const std::vector<int>& EdgeCut);

template<typename T>
bool ExportList(const std::string& Filename,
                const std::vector<T>& List)
{
    std::ofstream Stream;
    Stream.open(Filename, std::ios::out);
    if (!Stream.is_open())
        return false;

    for (const auto& v : List)
        Stream << v << '\n';

    Stream.close();
    return true;
}

template<typename T>
bool LoadList(const std::string& Filename,
              std::vector<T>& List)
{
    std::ifstream Stream;
    Stream.open(Filename, std::ios::in);
    if (!Stream.is_open())
        return false;

    std::string Line;
    while (!Stream.eof())
    {
        std::getline(Stream, Line);
        if (Line.empty())
            continue;
        std::istringstream ss(Line);
        T v;
        ss >> v;
        List.emplace_back(v);
    }

    Stream.close();
    return true;
}

} // namespace dfy
