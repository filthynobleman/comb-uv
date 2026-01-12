/**
 * @file        main.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#include <dfy/io.hpp>
#include <dfy/segmentation.hpp>
#include <dfy/sampler.hpp>
#include <iostream>

int main(int argc, const char* const argv[])
{
    dfy::ManifoldMesh M(argv[1]);
    std::cout << "Loaded mesh " << argv[1] << std::endl;
    dfy::Graph G = dfy::DualMeshToGraph(M);
    std::cout << "Converted to graph" << std::endl;
    dfy::Sampler Smpl(G);
    while (Smpl.NumSamples() < 10)
        Smpl.AddSample();
    std::cout << "Voronoi computed" << std::endl;
    std::vector<int> VParts(Smpl.GetPartitions().data(), 
                            Smpl.GetPartitions().data() + M.NumTriangles());
    dfy::ExportList("./voronoi.txt", VParts);
    std::cout << "Voronoi exported" << std::endl;
    dfy::Segmentation Seg(M, G, Smpl.GetPartitions());
    std::cout << "Segmentation built" << std::endl;
    Seg.MakeAllDisks(5);
    std::cout << "Disks computed" << std::endl;
    std::vector<int> DParts(Seg.GetTriParts().data(), 
                            Seg.GetTriParts().data() + M.NumTriangles());
    dfy::ExportList("./disks.txt", DParts);
    std::cout << "Disks exported" << std::endl;

    return 0;
}
