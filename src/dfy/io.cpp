/**
 * @file        io.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#include <dfy/io.hpp>

#include <igl/write_triangle_mesh.h>


bool dfy::ExportMesh(const std::string &Filename, 
                     const dfy::Mesh &M)
{
    return igl::write_triangle_mesh(Filename, M.Vertices(), M.Triangles());
}


bool dfy::ExportGraph(const std::string &Filename, 
                      const dfy::Graph &G, 
                      const Eigen::MatrixXd &V)
{
    std::ofstream Stream;
    Stream.open(Filename, std::ios::out);
    if (!Stream.is_open())
        return false;

    Stream << "o " << Filename.substr(0, Filename.rfind('.')) << '\n';

    for (int i = 0; i < G.NumNodes(); ++i)
        Stream << "v " << V(i, 0) << ' ' << V(i, 1) << ' ' << V(i, 2) << '\n';
    
    for (int i = 0; i < G.NumNodes(); ++i)
    {
        int nadj = G.NumAdjacents(i);
        for (int j = 0; j < nadj; ++j)
        {
            int adj = G.GetAdjacent(i, j).first;
            if (adj < i)
                continue;
            Stream << "l " << (i + 1) << ' ' << (adj + 1) << '\n';
        }
    }

    Stream.close();
    return true;
}