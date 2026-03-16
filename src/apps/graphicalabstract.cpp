/**
 * @file        uvpatches.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-19
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


enum GraphMetric
{
    EUCLIDEAN,
    GEODESIC,
    ANGULAR
};


struct 
{
    std::string InputFile;
    std::string OutputFile;
    int StartSamples;
    int SubSamples;
    double Threshold;
    GraphMetric Metric;
    bool SmoothNormals;

    int Verbosity;
} CLIArgs;


void ParseArgs(int argc, const char* const argv[]);

void StartTimer();
double StopTimer();
void StartGlobalTimer();
double StopGlobalTimer();


int main(int argc, const char* const argv[])
{
    ParseArgs(argc, argv);

    std::filesystem::path OPath(CLIArgs.OutputFile);

    std::filesystem::path VoronoiPath = OPath.filename().replace_extension("");
    VoronoiPath = std::filesystem::path(VoronoiPath.string() + "-voronoi.obj");
    VoronoiPath = OPath.parent_path() / VoronoiPath;

    std::filesystem::path DiskPath = OPath.filename().replace_extension("");
    DiskPath = std::filesystem::path(DiskPath.string() + "-disks.obj");
    DiskPath = OPath.parent_path() / DiskPath;

    std::filesystem::path CutPath = OPath.filename().replace_extension("");
    CutPath = std::filesystem::path(CutPath.string() + "-cut.obj");
    CutPath = OPath.parent_path() / CutPath;

    // Load mesh
    StartTimer();
    dfy::ManifoldMesh M(CLIArgs.InputFile);
    M.FixFaceOrientation();
    std::cout << "Loaded mesh " << CLIArgs.InputFile;
    std::cout << " in " << StopTimer() << " seconds.";
    std::cout << std::endl;

    // Start global timer
    StartGlobalTimer();

    // Compute dual graph
    StartTimer();
    dfy::DM2GDist GMetric = dfy::DualAngularDistance;
    switch (CLIArgs.Metric)
    {
    case GraphMetric::EUCLIDEAN:
        GMetric = dfy::DualEuclideanDistance;
        break;

    case GraphMetric::GEODESIC:
        GMetric = dfy::GeodesicDistance;
        break;
    
    default:
        break;
    }
    dfy::Graph G = dfy::DualMeshToGraph(M, GMetric);
    std::cout << "Dual graph computed in ";
    std::cout << StopTimer() << " seconds.";
    std::cout << std::endl;

    // Voronoi decomposition
    StartTimer();
    dfy::Sampler Smpl(G);
    CLIArgs.StartSamples = std::max(CLIArgs.StartSamples, 2 * (M.Genus() + 1));
    Smpl.AddSamples(std::min(CLIArgs.StartSamples, G.NumNodes()) - 1);
    std::cout << "Voronoi decomposition with ";
    std::cout << Smpl.NumSamples();
    std::cout << " samples computed in ";
    std::cout << StopTimer() << " seconds.";
    std::cout << std::endl;

    dfy::ExportMesh(VoronoiPath.string(), M.Separate(Smpl.GetPartitions()));

    // Disk decomposition
    StartTimer();
    dfy::Segmentation Seg(M, G, Smpl.GetPartitions());
    Seg.MakeAllDisks(CLIArgs.SubSamples);
    std::cout << "Disk decomposition into ";
    std::cout << Seg.NumRegions() << " regions ";
    std::cout << "computed in " << StopTimer() << " seconds.";
    std::cout << std::endl;

    dfy::ExportMesh(DiskPath.string(), M.Separate(Seg.GetTriParts()));

    // Compute cut
    StartTimer();
    std::vector<int> ECut;
    Seg.CutToDisk(ECut);
    dfy::Mesh CM = M.CutEdges(ECut);
    std::cout << "Cut along " << ECut.size();
    std::cout << " edges computed in ";
    std::cout << StopTimer() << " seconds.";
    std::cout << std::endl;

    dfy::ExportMesh(CutPath.string(), CM);


    // Merge disks
    StartTimer();
    if (CLIArgs.Metric != GraphMetric::ANGULAR)
        Seg.MergeRegions(dfy::MinPerimeterScore, CLIArgs.Threshold);
    else
        Seg.MergeRegions(dfy::MaxAvgCurvature, CLIArgs.Threshold);
    std::cout << "Merged into " << Seg.NumRegions();
    std::cout << " regions in ";
    std::cout << StopTimer() << " seconds.";
    std::cout << std::endl;

    dfy::ExportMesh(CLIArgs.OutputFile, M.Separate(Seg.GetTriParts()));

    double Runtime = StopGlobalTimer();
    std::cout << "Total runtime: ";
    std::cout << Runtime << " seconds.";
    std::cout << std::endl;

    return EXIT_SUCCESS;
}


void ParseArgs(int argc, const char *const argv[])
{
    CLIArgs.InputFile = "../samples/fertility.obj";
    CLIArgs.OutputFile = "./fertility-vororefine.obj";
    CLIArgs.Metric = GraphMetric::GEODESIC;
    CLIArgs.StartSamples = 5;
    CLIArgs.SubSamples = 5;
    CLIArgs.Threshold = 0.76;
    CLIArgs.SmoothNormals = true;
    CLIArgs.Verbosity = 2;
}


std::chrono::system_clock::time_point m_Start;
std::chrono::system_clock::time_point m_GlobalStart;
void StartTimer()
{
    m_Start = std::chrono::system_clock::now();
}
void StartGlobalTimer()
{
    m_GlobalStart = std::chrono::system_clock::now();
}

double StopTimer()
{
    auto End = std::chrono::system_clock::now();
    auto ETA = End - m_Start;
    size_t ETAms = std::chrono::duration_cast<std::chrono::milliseconds>(ETA).count();
    return 1.0e-3 * ETAms;
}
double StopGlobalTimer()
{
    auto End = std::chrono::system_clock::now();
    auto ETA = End - m_GlobalStart;
    size_t ETAms = std::chrono::duration_cast<std::chrono::milliseconds>(ETA).count();
    return 1.0e-3 * ETAms;
}