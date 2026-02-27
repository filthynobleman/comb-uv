/**
 * @file        uvpatches_experiment.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-02-25
 */
#include <dfy/manifoldmesh.hpp>
#include <dfy/sampler.hpp>
#include <dfy/segmentation.hpp>
#include <dfy/arap.hpp>
#include <dfy/harmonic.hpp>
#include <dfy/conformal.hpp>
#include <dfy/io.hpp>
#include <dfy/utils.hpp>

#include <iostream>
#include <chrono>
#include <filesystem>

#define _USE_MATH_DEFINES
#include <math.h>



int main(int argc, const char* const argv[])
{
    dfy::ManifoldMesh M(argv[1]);
    M.FixFaceOrientation();

    dfy::Graph G = dfy::DualMeshToGraph(M, dfy::DualAngularDistance);
    dfy::Sampler Smpl(G);
    Smpl.AddSamples(std::max(2 * (M.Genus() + 1), 4));
    dfy::Segmentation Seg(M, G, Smpl.GetPartitions());
    Seg.MakeAllDisks(5);
    Seg.MergeRegions(dfy::MaxAvgCurvature, std::stoi(argv[2]));

    for (int a = 1; a <= 3; ++a)
    {
        dfy::UVMapAlgorithm Algo = (dfy::UVMapAlgorithm)a;

        Eigen::MatrixXd UV;
        Eigen::MatrixXi FUV;
        FUV.setConstant(M.NumTriangles(), 3, -1);
        for (int i = 0; i < Seg.NumRegions(); ++i)
        {
            std::vector<int> FID;
            Seg.GetRegion(i).Faces(FID);
            dfy::Mesh SM = M.SubMesh(FID, dfy::SubMeshAccessType::BY_TRIANGLES);
            Eigen::MatrixXd LUV;
            if (Algo == dfy::UVMapAlgorithm::TUTTE)
            {
                dfy::TutteEmbedding Emb(SM);
                Emb.MapBoundary(dfy::BoundaryMap::SQUARE);
                if (!Emb.Compute())
                {
                    std::cerr << "TUTTE embedding failed." << std::endl;
                    continue;
                }
                LUV = Emb.UV();
            }
            else if (Algo == dfy::UVMapAlgorithm::HARMONIC)
            {
                dfy::HarmonicEmbedding Emb(SM);
                Emb.MapBoundary(dfy::BoundaryMap::CIRCLE);
                if (!Emb.Compute())
                {
                    std::cerr << "HARMONIC embedding failed." << std::endl;
                    continue;
                }
                LUV = Emb.UV();
            }
            else if (Algo == dfy::UVMapAlgorithm::CONFORMAL)
            {
                dfy::ConformalEmbedding Emb(SM);
                if (!Emb.Compute())
                {
                    std::cerr << "CONFORMAL embedding failed." << std::endl;
                    continue;
                }
                LUV = Emb.UV();
            }
            else if (Algo == dfy::UVMapAlgorithm::ARAP)
            {
                dfy::ARAPEmbedding Emb(SM);
                if (!Emb.Compute())
                {
                    std::cerr << "ARAP embedding failed." << std::endl;
                    continue;
                }
                LUV = Emb.UV();
            }
            int NUVs = UV.rows();
            UV.conservativeResize(NUVs + LUV.rows(), 2);
            UV(Eigen::seq(NUVs, UV.rows() - 1), Eigen::all) = LUV;
            FUV(FID, Eigen::all) = SM.Triangles().array() + NUVs;
        }

        dfy::ExportMesh(argv[2 + a], M, UV, FUV, true);
    }

    return EXIT_SUCCESS;
}