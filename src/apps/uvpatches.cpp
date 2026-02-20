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

    // Load mesh
    StartTimer();
    dfy::ManifoldMesh M(CLIArgs.InputFile);
    M.FixFaceOrientation();
    if (CLIArgs.Verbosity > 1)
    {
        std::cout << "Loaded mesh " << CLIArgs.InputFile;
        std::cout << " in " << StopTimer() << " seconds.";
        std::cout << std::endl;
    }

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
    if (CLIArgs.Verbosity > 1)
    {
        std::cout << "Dual graph computed in ";
        std::cout << StopTimer() << " seconds.";
        std::cout << std::endl;
    }

    // Voronoi decomposition
    StartTimer();
    dfy::Sampler Smpl(G);
    CLIArgs.StartSamples = std::max(CLIArgs.StartSamples, 4 * M.Genus());
    Smpl.AddSamples(std::min(CLIArgs.StartSamples, G.NumNodes()) - 1);
    if (CLIArgs.Verbosity > 1)
    {
        std::cout << "Voronoi decomposition with ";
        std::cout << Smpl.NumSamples();
        std::cout << " samples computed in ";
        std::cout << StopTimer() << " seconds.";
        std::cout << std::endl;
    } 

    // Disk decomposition
    StartTimer();
    dfy::Segmentation Seg(M, G, Smpl.GetPartitions());
    Seg.MakeAllDisks(CLIArgs.SubSamples);
    if (CLIArgs.Verbosity > 1)
    {
        std::cout << "Disk decomposition computed in ";
        std::cout << StopTimer() << " seconds.";
        std::cout << std::endl;
    }

    // Merge disks
    StartTimer();
    if (CLIArgs.Metric != GraphMetric::ANGULAR)
        Seg.MergeRegions(dfy::MinPerimeterScore, CLIArgs.Threshold);
    else
        Seg.MergeRegions(dfy::MaxAvgDihedralAngle, CLIArgs.Threshold);
    if (CLIArgs.Verbosity > 1)
    {
        std::cout << "Merged regions in ";
        std::cout << StopTimer() << " seconds.";
        std::cout << std::endl;
    }

    // Unwrap each region
    StartTimer();
    Eigen::MatrixXd UV;
    Eigen::MatrixXi FUV;
    FUV.setConstant(M.NumTriangles(), 3, -1);
    for (int i = 0; i < Seg.NumRegions(); ++i)
    {
        std::vector<int> FID;
        Seg.GetRegion(i).Faces(FID);
        dfy::Mesh SM = M.SubMesh(FID, dfy::SubMeshAccessType::BY_TRIANGLES);
        Eigen::MatrixXd LUV;
        if (CLIArgs.Algorithm == dfy::UVMapAlgorithm::TUTTE)
        {
            dfy::TutteEmbedding Emb(SM);
            Emb.MapBoundary(dfy::BoundaryMap::SQUARE);
            Emb.Compute();
            LUV = Emb.UV();
        }
        else if (CLIArgs.Algorithm == dfy::UVMapAlgorithm::HARMONIC)
        {
            dfy::HarmonicEmbedding Emb(SM);
            Emb.MapBoundary(dfy::BoundaryMap::CIRCLE);
            Emb.Compute();
            LUV = Emb.UV();
        }
        else if (CLIArgs.Algorithm == dfy::UVMapAlgorithm::CONFORMAL)
        {
            dfy::ConformalEmbedding Emb(SM);
            Emb.Compute();
            LUV = Emb.UV();
        }
        else if (CLIArgs.Algorithm == dfy::UVMapAlgorithm::ARAP)
        {
            dfy::ARAPEmbedding Emb(SM);
            Emb.Compute();
            LUV = Emb.UV();
        }
        if (CLIArgs.Packing)
        {
            int GSize = std::ceil(std::sqrt(Seg.NumRegions()));
            LUV = (LUV.array() * 0.975 + 0.0125) / GSize;
            LUV.col(0) = LUV.col(0).array() + (i % GSize) / ((double) GSize);
            LUV.col(1) = LUV.col(1).array() + (i / GSize) / ((double) GSize);
        }
        int NUVs = UV.rows();
        UV.conservativeResize(NUVs + LUV.rows(), 2);
        UV(Eigen::seq(NUVs, UV.rows() - 1), Eigen::all) = LUV;
        FUV(FID, Eigen::all) = SM.Triangles().array() + NUVs;
    }
    if (CLIArgs.Verbosity > 1)
    {
        std::cout << "Computed embedding of ";
        std::cout << Seg.NumRegions() << " patches in ";
        std::cout << StopTimer() << " seconds.";
        std::cout << std::endl;
    }

    double Runtime = StopGlobalTimer();

    StartTimer();
    dfy::ExportMesh(CLIArgs.OutputFile, M, UV, FUV, CLIArgs.SmoothNormals);
    if (CLIArgs.Verbosity > 1)
    {
        std::cout << "Mesh exported to " << CLIArgs.OutputFile << " in ";
        std::cout << StopTimer() << " seconds.";
        std::cout << std::endl;
    }

    
    if (CLIArgs.Verbosity > 0)
    {
        std::cout << "Total runtime: ";
        std::cout << Runtime << " seconds.";
        std::cout << std::endl;
    }

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