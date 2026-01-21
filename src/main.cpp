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
#include <dfy/arap.hpp>
#include <dfy/gimage.hpp>
#include <dfy/imgsampler.hpp>
#include <iostream>

int main(int argc, const char* const argv[])
{
    dfy::ManifoldMesh M(argv[1]);
    std::cout << "Loaded mesh " << argv[1] << std::endl;
    dfy::Graph G = dfy::DualMeshToGraph(M, dfy::DualAngularDistance);
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
    std::cout << "Disks computed (" << Seg.NumRegions() << ")" << std::endl;
    std::vector<int> DParts(Seg.GetTriParts().data(), 
                            Seg.GetTriParts().data() + M.NumTriangles());
    dfy::ExportList("./disks.txt", DParts);
    std::cout << "Disks exported" << std::endl;
    std::vector<int> CutEdges;
    Eigen::VectorXd EWeights;
    EWeights.resize(M.NumEdges());
    for (int e = 0; e < M.NumEdges(); ++e)
    {
        int i = M.EdgeTriAdj()(e, 0);
        int j = M.EdgeTriAdj()(e, 1);
        // int i = M.Edges()(e, 0);
        // int j = M.Edges()(e, 1);
        if (i == -1 || j == -1)
            EWeights[e] = 0;
        else
            EWeights[e] = dfy::DualAngularDistance(M, i, j);
    }
    EWeights = EWeights.maxCoeff() - EWeights.array();
    Seg.CutToDisk(CutEdges, EWeights);
    std::cout << "Cut to disk computed" << std::endl;
    dfy::ExportPointCloud("./edges.obj", M.EdgeCenters()(CutEdges, Eigen::all).eval());
    dfy::Mesh DM = M.CutEdges(CutEdges);
    std::cout << "Cut executed" << std::endl;
    dfy::ExportMesh("./cut.obj", DM);
    std::cout << "Cut exported" << std::endl;
    dfy::TutteEmbedding Emb(DM);
    std::cout << "Boundary length: " << Emb.BoundaryLength() << std::endl;
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
    dfy::ImageSampler ImgSmpl(Img);
    Eigen::MatrixXd QVerts;
    int QSize = 256;
    QVerts.resize(QSize * QSize, 3);
    for (int j = 0; j < QSize; ++j)
    {
        double v = j / (QSize - 1.0);
        for (int i = 0; i < QSize; ++i)
        {
            double u = i / (QSize - 1.0);
            int RowIdx = j * QSize + i;
            QVerts.row(RowIdx) = ImgSmpl.Sample(u, v, dfy::ImgInterpType::CUBIC);
            // QVerts.row(RowIdx) = Img.GetPixel(i, j);
        }
    }
    std::cout << "Sampled image" << std::endl;
    dfy::QuadMesh QM(QVerts, QSize, QSize);
    std::cout << "Created quad mesh" << std::endl;
    dfy::ExportQuadMesh("./qmesh.obj", QM);
    std::cout << "Exported quad mesh" << std::endl;

    return EXIT_SUCCESS;
}
