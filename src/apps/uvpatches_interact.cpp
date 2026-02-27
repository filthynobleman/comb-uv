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

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

// #include <imgui/imgui.h>

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
    dfy::UVMapAlgorithm Algorithm;
    GraphMetric Metric;
    bool Packing;
    bool SmoothNormals;

    int Verbosity;
} CLIArgs;


void ParseArgs(int argc, const char* const argv[]);
void Usage(const std::string& argv0);

void StartTimer();
double StopTimer();
void StartGlobalTimer();
double StopGlobalTimer();


int main(int argc, const char* const argv[])
{
    ParseArgs(argc, argv);

    polyscope::init();

    dfy::ManifoldMesh M(CLIArgs.InputFile);

    std::string MeshName = std::filesystem::path(CLIArgs.InputFile).filename().replace_extension("").string();
    auto psMesh = polyscope::registerSurfaceMesh(MeshName, M.Vertices(), M.Triangles());

    int NSamples = 2 * (M.Genus() + 1);
    dfy::Graph G = dfy::DualMeshToGraph(M, dfy::DualAngularDistance);
    dfy::Sampler Smpl(G);
    Smpl.AddSamples(NSamples - 1);
    dfy::Segmentation Seg(M, G, Smpl.GetPartitions());
    std::vector<Eigen::Vector3d> VoronoiColors;
    VoronoiColors.resize(M.NumTriangles());
    for (int i = 0; i < Seg.NumRegions(); ++i)
    {
        std::vector<int> FID;
        Seg.GetRegion(i).Faces(FID);
        Eigen::Vector3d CurCol{ polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit() };
        for (int f : FID)
            VoronoiColors[f] = CurCol;
    }
    Seg.MakeAllDisks(NSamples);
    std::vector<Eigen::Vector3d> DisksColor;
    DisksColor.resize(M.NumTriangles());
    for (int i = 0; i < Seg.NumRegions(); ++i)
    {
        std::vector<int> FID;
        Seg.GetRegion(i).Faces(FID);
        Eigen::Vector3d CurCol{ polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit() };
        for (int f : FID)
            DisksColor[f] = CurCol;
    }
    std::vector<Eigen::Vector3d> MergeColor;
    MergeColor = DisksColor;
    float Threshold = 1e-2;

    // MergeColor.resize(M.NumTriangles());
    // for (int i = 0; i < Seg.NumRegions(); ++i)
    // {
    //     std::vector<int> FID;
    //     Seg.GetRegion(i).Faces(FID);
    //     Eigen::Vector3d CurCol{ polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit() };
    //     for (int f : FID)
    //         MergeColor[f] = CurCol;
    // }

    psMesh->addFaceColorQuantity("Voronoi", VoronoiColors);
    psMesh->addFaceColorQuantity("Disks", DisksColor);
    auto MrgQnt = psMesh->addFaceColorQuantity("Merge", MergeColor);

    auto MergeCbk = [&]() {
        ImGui::InputFloat("Threshold", &Threshold, 0.0f, 0.0f, "%.3e");
        if (!ImGui::Button("Merge regions"))
            return;
        
        dfy::Segmentation SegMerge(Seg);

        SegMerge.MergeRegions(dfy::MaxAvgCurvature, Threshold);
        MergeColor.resize(M.NumTriangles());
        for (int i = 0; i < SegMerge.NumRegions(); ++i)
        {
            std::vector<int> FID;
            SegMerge.GetRegion(i).Faces(FID);
            Eigen::Vector3d CurCol{ polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit() };
            for (int f : FID)
                MergeColor[f] = CurCol;
        }

        MrgQnt->updateData(MergeColor);
    };

    polyscope::state::userCallback = MergeCbk;
    polyscope::show();

    return EXIT_SUCCESS;
}


void ParseArgs(int argc, const char *const argv[])
{
    CLIArgs.InputFile = "../samples/bunny.obj";
    CLIArgs.Algorithm = dfy::UVMapAlgorithm::TUTTE;
    CLIArgs.Metric = GraphMetric::GEODESIC;
    CLIArgs.StartSamples = 5;
    CLIArgs.SubSamples = 5;
    CLIArgs.Threshold = 0.76;
    CLIArgs.Packing = false;
    CLIArgs.SmoothNormals = false;
    CLIArgs.Verbosity = 1;

    for (int i = 1; i < argc; ++i)
    {
        if (argv[i][0] != '-')
        {
            CLIArgs.InputFile = argv[i];
            continue;
        }

        std::string Argvi = argv[i];
        // UV packing
        if (Argvi == "-p" || Argvi == "--packing")
        {
            CLIArgs.Packing = true;
            continue;
        }
        // Smooth normals in output
        if (Argvi == "-s" || Argvi == "--smooth")
        {
            CLIArgs.SmoothNormals = true;
            continue;
        }
        // Help
        if (Argvi == "-h" || Argvi == "--help")
        {
            Usage(argv[0]);
            exit(EXIT_SUCCESS);
        }

        // Other arguments require specifying a value
        if (i + 1 == argc)
        {
            std::cerr << "Option \"" << Argvi << "\" has no specified value.";
            std::cerr << std::endl;
            Usage(argv[0]);
            exit(EXIT_FAILURE);
        }
        // Starting samples for Voronoi
        if (Argvi == "-n" || Argvi == "--num-samples")
        {
            CLIArgs.StartSamples = std::stoi(argv[++i]);
            continue;
        }
        // Subsamples for disk decomposition
        if (Argvi == "-d" || Argvi == "--disk-samples")
        {
            CLIArgs.SubSamples = std::stoi(argv[++i]);
            continue;
        }
        // Merging threshold
        if (Argvi == "-t" || Argvi == "--threshold")
        {
            CLIArgs.Threshold = std::stod(argv[++i]);
            continue;
        }
        // Unwrapping algorithm
        if (Argvi == "-a" || Argvi == "--algorithm")
        {
            std::string Algo = argv[++i];
            std::transform(Algo.begin(), Algo.end(), Algo.begin(),
                           [](unsigned char c) { return std::tolower(c); });

            if (Algo == "tutte")
                CLIArgs.Algorithm = dfy::UVMapAlgorithm::TUTTE;
            else if (Algo == "harmonic")
                CLIArgs.Algorithm = dfy::UVMapAlgorithm::HARMONIC;
            else if (Algo == "arap")
                CLIArgs.Algorithm = dfy::UVMapAlgorithm::ARAP;
            else if (Algo == "conformal")
                CLIArgs.Algorithm = dfy::UVMapAlgorithm::CONFORMAL;
            else
            {
                std::cerr << "Specified unwrapping algorithm \"";
                std::cerr << Algo << "\" is not a valid option.";
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
                Usage(argv[0]);
            }
        }
        // Graph metric
        if (Argvi == "-m" || Argvi == "--metric")
        {
            std::string Metric = argv[++i];
            std::transform(Metric.begin(), Metric.end(), Metric.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (Metric == "angular")
                CLIArgs.Metric = GraphMetric::ANGULAR;
            else if (Metric == "geodesic")
                CLIArgs.Metric = GraphMetric::GEODESIC;
            else if (Metric == "euclidean")
                CLIArgs.Metric = GraphMetric::EUCLIDEAN;
            else
            {
                std::cerr << "Specified graph metric \"";
                std::cerr << Metric << "\" is not a valid option.";
                std::cerr << std::endl;
                Usage(argv[0]);
                exit(EXIT_FAILURE);
            }
        }
        // Output mesh
        if (Argvi == "-o" || Argvi == "--output")
        {
            CLIArgs.OutputFile = argv[++i];
            continue;
        }
        // Verbosity level
        if (Argvi == "-v" || Argvi == "--verbosity")
        {
            CLIArgs.Verbosity = std::stoi(argv[++i]);
            continue;
        }
    }

    // If no output mesh is given, use input mesh with "-uv" as suffix
    if (CLIArgs.OutputFile.empty())
    {
        std::filesystem::path InPath(CLIArgs.InputFile);
        InPath.replace_filename(InPath.stem().string() + "-uv.obj");
        CLIArgs.OutputFile = InPath.string();
    }

    // If extension is not obj, raise a warning
    std::filesystem::path OutPath(CLIArgs.OutputFile);
    std::string Ext = OutPath.extension().string();
    std::transform(Ext.begin(), Ext.end(), Ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    OutPath.replace_extension(".obj");
    if (Ext != ".obj")
    {
        std::cerr << "Output filename has no OBJ extension (";
        std::cerr << CLIArgs.OutputFile << "). ";
        std::cerr << "Changing name to " << OutPath.string();
        std::cerr << std::endl;
        CLIArgs.OutputFile = OutPath.string();
    }
}

void Usage(const std::string &argv0)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "  " << argv0 << " input_mesh [-n num_samples] [-d disk_samples] [-t threshold] [-o output_mesh] [-a algorithm] [-v verbosity] [-p] [-s]" << std::endl;
    std::cout << "  " << argv0 << " input_mesh [--num-samples num_samples] [--disk-samples disk_samples] [--threshold threshold] [--output output_mesh] [--algorithm algorithm] [--verbosity verbosity] [--packing] [--smooth]" << std::endl;
    std::cout << "  " << argv0 << " -h" << std::endl;
    std::cout << "  " << argv0 << " --help" << std::endl;
    std::cout << std::endl;
    std::cout << "    input_mesh is the path to a triangulated mesh." << std::endl;
    std::cout << "    num_samples is the number of samples for the initial Voronoi partitioning. Default is 5." << std::endl;
    std::cout << "    disk_samples is the number of sub-samples that must subdivide non topological disks. Default is 5." << std::endl;
    std::cout << "    threshold is the threshold parameter that determines if two regions can be merged. Default is 0.76" << std::endl;
    std::cout << "    algorithm is the UV unwrapping algorithm for each region. Acceptable values are \'tutte\', \'harmonic\', \'conformal\', \'arap\'. Default is \'harmonic\'." << std::endl;
    std::cout << "    verbosity is the verbosity level of the output, ranging from 0 (no output) to 2 (runtime of each step). Default is 1 (total runtime)." << std::endl;
    std::cout << "    -p|--packing option forces each UV island to not overlap with the others. By default, island are all rescaled to [0, 1]^2." << std::endl;
    std::cout << "    -s|--smooth option outputs a mesh with smoothed normals. By default, output mesh has constant normals over triangles." << std::endl;
    std::cout << "    -h|--help option prints this help message." << std::endl;
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