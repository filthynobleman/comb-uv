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
#include <dfy/tutte.hpp>
#include <dfy/gimage.hpp>
#include <iostream>

int main(int argc, const char* const argv[])
{
    dfy::ManifoldMesh M(argv[1]);
    std::cout << "Loaded mesh " << argv[1] << std::endl;
    dfy::Graph G = dfy::DualMeshToGraph(M, dfy::DualCurvatureDistance);
    std::cout << "Converted to graph (" << G.NumEdges() / 2 << " edges)" << std::endl;
    dfy::Sampler Smpl(G);
    Smpl.AddSamples(5);
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
    std::vector<int> CutEdges;
    Seg.CutToDisk(CutEdges);
    std::cout << "Cut to disk computed" << std::endl;
    dfy::ExportPointCloud("./edges.obj", M.EdgeCenters()(CutEdges, Eigen::all).eval());
    dfy::Mesh DM = M.CutEdges(CutEdges);
    std::cout << "Cut executed" << std::endl;
    dfy::ExportMesh("./cut.obj", DM);
    std::cout << "Cut exported" << std::endl;
    dfy::TutteEmbedding Emb(DM);
    Emb.MapBoundary(dfy::BoundaryMap::SQUARE);
    if (!Emb.Compute())
    {
        std::cerr << "Something went wrong" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Embedding computed" << std::endl;
    if (!dfy::ExportMesh("./unwrapped.obj", M, Emb.UV(), DM.Triangles()))
    {
        std::cerr << "Something went wrong" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Mesh exported" << std::endl;
    dfy::GImage Img(M, 4096);
    std::cout << "Initialized geometry image" << std::endl;
    Img.Compute(Emb.UV(), DM.Triangles());
    std::cout << "Computed geometry image" << std::endl;
    dfy::ExportGImage("./gimage.png", Img);
    std::cout << "Exported geometry image" << std::endl;

    return EXIT_SUCCESS;
}
