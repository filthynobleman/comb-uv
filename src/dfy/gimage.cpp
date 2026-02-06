/**
 * @file        gimage.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-14
 */
#include <dfy/gimage.hpp>

#include <dfy/utils.hpp>
#include <dfy/imgsampler.hpp>

#include <igl/barycentric_coordinates.h>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

#include <fstream>

dfy::GImage::GImage(const dfy::Mesh &M, 
                    int Width, 
                    int Height)
    : m_Mesh(M), m_Width(Width), m_Height(Height)
{
    m_Data = (Eigen::Vector3d*)std::calloc(Width * Height, sizeof(Eigen::Vector3d));
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot create GImage.");
}

dfy::GImage::GImage(const dfy::Mesh &M, int Size)
    : dfy::GImage(M, Size, Size) { }

dfy::GImage::GImage(const dfy::GImage &GI)
    : m_Mesh(GI.m_Mesh)
{
    m_Width = GI.m_Width;
    m_Height = GI.m_Height;
    size_t BufSize = m_Width * m_Height * sizeof(Eigen::Vector3d);
    m_Data = (Eigen::Vector3d*)std::malloc(BufSize);
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot create GImage.");
    m_Data = (Eigen::Vector3d*)std::memcpy(m_Data, GI.m_Data, BufSize);
    if (m_Data == nullptr)
        throw std::runtime_error("Cannot create GImage.");
}

dfy::GImage::GImage(dfy::GImage &&GI)
    : m_Mesh(GI.m_Mesh)
{
    m_Width = GI.m_Width;
    m_Height = GI.m_Height;
    m_Data = GI.m_Data;
    GI.m_Data = nullptr;
}

dfy::GImage::~GImage()
{
    delete m_Data;
}

int dfy::GImage::GetWidth() const { return m_Width; }
int dfy::GImage::GetHeight() const { return m_Height; }

const Eigen::Vector3d &dfy::GImage::GetPixel(int i, int j) const
{
    return m_Data[j * GetWidth() + i];
}

void dfy::GImage::GetPixel(int i, int j, 
                           Eigen::Vector3d &rgb) const
{
    rgb = GetPixel(i, j);
}

void dfy::GImage::GetPixel(int i, int j, 
                           unsigned char &r, unsigned char &g, unsigned char &b) const
{
    const Eigen::Vector3d& rgb = GetPixel(i, j);
    r = std::round(std::min(1.0, std::max(0.0, rgb.x())) * 255.0);
    g = std::round(std::min(1.0, std::max(0.0, rgb.y())) * 255.0);
    b = std::round(std::min(1.0, std::max(0.0, rgb.z())) * 255.0);
}

#include <iostream>
#include <random>

void dfy::GImage::Compute(const Eigen::MatrixXd &UV, 
                          const Eigen::MatrixXi &TUV)
{
    // Normalize mesh coordinates
    Eigen::MatrixXd RGB = m_Mesh.Vertices();
    for (int i = 0; i < 3; ++i)
        RGB.col(i) = RGB.col(i).array() - RGB.col(i).minCoeff();
    RGB /= RGB.maxCoeff();
    // Check each pixel has been set
    std::vector<bool> PixSet;
    PixSet.resize(GetWidth() * GetHeight(), false);
    for (int t = 0; t < TUV.rows(); ++t)
    {
        // Get triangle
        Eigen::Matrix<double, 3, 2> Tri = UV(TUV(t, Eigen::all).transpose(), Eigen::all);
        // Get bounding box of triangle
        Eigen::Vector2d dBBL = Tri.colwise().minCoeff().cwiseMax(0.0);
        Eigen::Vector2d dBTR = Tri.colwise().maxCoeff().cwiseMin(1.0);
        Eigen::Vector2i BBL, BTR;
        BBL[0] = std::floor(dBBL[0] * GetWidth());
        BBL[1] = std::floor(dBBL[1] * GetHeight());
        BTR[0] = std::ceil(dBTR[0] * GetWidth());
        BTR[1] = std::ceil(dBTR[1] * GetHeight());
        
        // Iterate over the pixels inside the bounding box
        for (int j = BBL.y(); j < BTR.y(); ++j)
        {
            for (int i = BBL.x(); i < BTR.x(); ++i)
            {
                // Get the pixel coordinates
                Eigen::RowVector2d ij{ (i + 1e-6) / (GetWidth() - (1 - 2e-6)), (j + 1e-6) / (GetHeight() - (1 - 2e-6)) };
                Eigen::RowVector3d Lambda;
                igl::barycentric_coordinates(ij, Tri.row(0), Tri.row(1), Tri.row(2), Lambda);
                if ((Lambda.array() > -1e-6).all())
                {
                    Eigen::Matrix3d T3D = RGB(m_Mesh.Triangles()(t, Eigen::all), Eigen::all);
                    m_Data[j * GetWidth() + i] = Lambda * T3D;
                    PixSet[j * GetWidth() + i] = true;
                }
            }
        }
    }

    // Fill unset pixels with closest triangle
    for (int j = 0; j < GetHeight(); ++j)
    {
        for (int i = 0; i < GetWidth(); ++i)
        {
            if (PixSet[j * GetWidth() + i])
                continue;

            // Get the pixel coordinates
            Eigen::RowVector2d ij{ (i + 1e-6) / (GetWidth() - (1 - 2e-6)), 
                                   (j + 1e-6) / (GetHeight() - (1 - 2e-6)) };
            
            // Find closest triangle
            int Closest = -1;
            double Error = std::numeric_limits<double>::infinity();
            for (int t = 0; t < TUV.rows(); ++t)
            {
                // Get triangle
                Eigen::Matrix<double, 3, 2> Tri = UV(TUV(t, Eigen::all).transpose(), Eigen::all);
                // Get barycentric coords
                Eigen::RowVector3d Lambda;
                igl::barycentric_coordinates(ij, Tri.row(0), Tri.row(1), Tri.row(2), Lambda);
                Lambda = Lambda.cwiseMax(0);
                Lambda /= Lambda.sum();
                double err = (ij - Lambda * Tri).eval().squaredNorm();
                if (err < Error)
                {
                    Error = err;
                    Closest = t;
                }
            }

            // Apply closest triangle
            int t = Closest;
            // Get triangle
            Eigen::Matrix<double, 3, 2> Tri = UV(TUV(t, Eigen::all).transpose(), Eigen::all);
            // Get barycentric coords
            Eigen::RowVector3d Lambda;
            igl::barycentric_coordinates(ij, Tri.row(0), Tri.row(1), Tri.row(2), Lambda);
            Lambda = Lambda.cwiseMax(0);
            Lambda /= Lambda.sum();
            Eigen::Matrix3d T3D = RGB(m_Mesh.Triangles()(t, Eigen::all), Eigen::all);
            m_Data[j * GetWidth() + i] = Lambda * T3D;
        }
    }
}

dfy::QuadMesh dfy::GImage::AsQuadMesh(const dfy::Mesh& M, int Res, bool Weld) const
{
    return AsQuadMesh(M, Res, Res, Weld);
}

// #include <igl/point_mesh_squared_distance.h>
// #include <igl/slim.h>

dfy::QuadMesh dfy::GImage::AsQuadMesh(const dfy::Mesh& M, int URes, int VRes, bool Weld) const
{
    // Initialize an image sampler
    dfy::ImageSampler ImgSampler(*this);

    // Initialize grid
    Eigen::MatrixXd U(VRes, URes);
    Eigen::MatrixXd V(VRes, URes);
    Eigen::MatrixXd UV(URes * VRes, 2);
    std::vector<int> BLoop;
    BLoop.reserve(2 * URes + 2 * VRes);
    std::vector<Eigen::Vector3i> Tris;
    Tris.reserve(2 * URes * VRes);
    for (int j = 0; j < VRes; ++j)
    {
        double v = j / (VRes - 1.0);
        for (int i = 0; i < URes; ++i)
        {
            double u = i / (URes - 1.0);
            U(j, i) = u;
            V(j, i) = v;
            UV.row(j * URes + i) = Eigen::Vector2d{ u, v };
            if (i == 0 || j == 0 || i == (URes - 1) || j == (VRes - 1))
                BLoop.emplace_back(j * URes + i);
            if (i < URes - 1 && j < VRes - 1)
            {
                int v1 = j * URes + i;
                int v2 = v1 + URes;
                int v3 = v2 + 1;
                int v4 = v1 + 1;
                Tris.emplace_back(v1, v2, v3);
                Tris.emplace_back(v1, v3, v4);
            }
        }
    }

    Eigen::VectorXi BL = Eigen::VectorXi::Map(BLoop.data(), BLoop.size());
    Eigen::MatrixXi Faces;
    Faces.resize(Tris.size(), 3);
    for (int i = 0; i < Tris.size(); ++i)
        Faces.row(i) = Tris[i];

    // Compute surface area
    double MeshArea = M.FaceAreas().sum();
    
    // Optimize
    const int MaxIter = 0;
    double Err = std::numeric_limits<double>::infinity();
    Eigen::MatrixXd Pos3D(URes * VRes, 3);
    Eigen::MatrixXd NewPos3D(URes * VRes, 3);
    #pragma omp parallel for
    for (int j = 0; j < VRes; ++j)
    {
        for (int i = 0; i < URes; ++i)
        {
            int Idx = j * URes + i;
            Pos3D.row(Idx) = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
        }
    }
    Eigen::MatrixXd BC = UV(BL, Eigen::all);
    NewPos3D = Pos3D;
    for (int Iter = 0; Iter < MaxIter && Err > 1e-4; ++Iter)
    {
        // // Update positions of interior vertices
        // igl::SLIMData SData;
        // igl::slim_precompute(Pos3D, Faces, UV, SData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, BL, BC, 0);
        // SData.mesh_area = MeshArea;
        // UV = igl::slim_solve(SData, 1);
        // #pragma omp parallel for
        // for (int j = 1; j < VRes - 1; ++j)
        // {
        //     for (int i = 1; i < URes - 1; ++i)
        //     {
        //         int Idx = j * URes + i;
        //         Pos3D.row(Idx) = ImgSampler.Sample(UV.row(Idx), dfy::ImgInterpType::CUBIC);
        //     }
        // }
        #pragma omp parallel for
        for (int j = 1; j < VRes - 1; ++j)
        {
            for (int i = 1; i < URes - 1; ++i)
            {
                int Idx = j * URes + i;
                NewPos3D.row(Idx).setZero();
                NewPos3D.row(Idx) += 0.25 * Pos3D.row(Idx - 1);
                NewPos3D.row(Idx) += 0.25 * Pos3D.row(Idx + 1);
                NewPos3D.row(Idx) += 0.25 * Pos3D.row(Idx - URes);
                NewPos3D.row(Idx) += 0.25 * Pos3D.row(Idx + URes);
            }
        }

        // // Project back onto the surface
        // Eigen::VectorXi I;
        // Eigen::VectorXd D;
        // Eigen::MatrixXd Tmp;
        // igl::point_mesh_squared_distance(NewPos3D, M.Vertices(), M.Triangles(), D, I, Tmp);
        const double alpha = 1e-4;
        #pragma omp parallel for
        for (int j = 1; j < VRes - 1; ++j)
        {
            for (int i = 1; i < URes - 1; ++i)
            {
                int Idx = j * URes + i;
                Eigen::Vector2d uv{ U(j, i), V(j, i) };
                Eigen::Vector3d p = ImgSampler.Sample(uv, dfy::ImgInterpType::CUBIC);
                Eigen::Matrix<double, 3, 2> J;
                // Projection by Newton iteration
                for (int k = 0; k < 10; ++k)
                {
                    J = ImgSampler.Jacobian(uv, dfy::ImgInterpType::CUBIC);
                    Eigen::Vector2d duv = (J.transpose() * J).eval().inverse() * (J.transpose() * (p - NewPos3D.row(Idx).transpose())).eval();
                    duv.normalize();
                    uv -= alpha * duv;
                    p = ImgSampler.Sample(uv, dfy::ImgInterpType::CUBIC);
                }
                U(j, i) = uv.x(); // 0.99 * U(j, i) + 0.01 * uv.x();
                V(j, i) = uv.y(); // 0.99 * V(j, i) + 0.01 * uv.y();
                Pos3D.row(Idx) = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
            }
        }

        // Compute residual
        double NewErr = 0.0;
        for (int j = 0; j < VRes; ++j)
        {
            for (int i = 0; i < URes; ++i)
            {
                int Idx = j * URes + 1;
                Eigen::Vector3d Pij = Pos3D.row(Idx);
                if (i < URes - 1)
                {
                    Eigen::Vector3d Po = Pos3D.row(Idx + 1);
                    NewErr += (Pij - Po).squaredNorm();
                }
                if (j < VRes - 1)
                {
                    Eigen::Vector3d Po = Pos3D.row(Idx + URes);
                    NewErr += (Pij - Po).squaredNorm();
                }
            }
        }
        if (std::abs(NewErr - Err) < 1e-6)
        {
            std::cout << "Converged to " << NewErr << " after " << Iter + 1 << " iterations." << std::endl;
            break;
        }
        Err = NewErr;
        std::cout << Err << std::endl;
    }
    if (Err > 1e-4)
        std::cerr << "WARNING: Optimization did not converge." << std::endl;
    std::cout << Err << std::endl;

    // // Create the quad mesh
    // Eigen::MatrixXd Verts(URes * VRes, 3);
    // int CurIdx = 0;
    // for (int j = 0; j < VRes; ++j)
    // {
    //     for (int i = 0; i < URes; ++i)
    //         Verts.row(CurIdx++) = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
    // }

    if (Weld)
        throw std::runtime_error("Boundary welding is not yet implemented.");

    return dfy::QuadMesh(Pos3D, URes, VRes);
}

// dfy::QuadMesh dfy::GImage::AsQuadMesh(int Res, bool Weld) const
// {
//     return AsQuadMesh(Res, Res, Weld);
// }

// dfy::QuadMesh dfy::GImage::AsQuadMesh(int URes, int VRes, bool Weld) const
// {
//     // Initialize an image sampler
//     dfy::ImageSampler ImgSampler(*this);

//     // Initialize grid
//     Eigen::MatrixXd U(VRes, URes);
//     Eigen::MatrixXd V(VRes, URes);
//     for (int j = 0; j < VRes; ++j)
//     {
//         double v = j / (VRes - 1.0);
//         for (int i = 0; i < URes; ++i)
//         {
//             double u = i / (URes - 1.0);
//             U(j, i) = u;
//             V(j, i) = v;
//         }
//     }

    
    
//     // Optimize
//     const int MaxIter = 100;
//     double Err = std::numeric_limits<double>::infinity();
//     std::vector<Eigen::Vector3d> Pos3D;
//     std::vector<Eigen::Vector3d> NewPos3D;
//     Pos3D.resize(URes * VRes);
//     #pragma omp parallel for
//     for (int j = 0; j < VRes; ++j)
//     {
//         for (int i = 0; i < URes; ++i)
//         {
//             int Idx = j * URes + i;
//             Pos3D[Idx] = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//         }
//     }
//     NewPos3D = Pos3D;
//     for (int Iter = 0; Iter < MaxIter && Err > 1e-4; ++Iter)
//     {
//         // Update positions of interior vertices
//         #pragma omp parallel for
//         for (int j = 1; j < VRes - 1; ++j)
//         {
//             for (int i = 1; i < URes - 1; ++i)
//             {
//                 int Idx = j * URes + i;
//                 NewPos3D[Idx].setZero();
//                 NewPos3D[Idx] += 0.25 * Pos3D[Idx - 1];
//                 NewPos3D[Idx] += 0.25 * Pos3D[Idx + 1];
//                 NewPos3D[Idx] += 0.25 * Pos3D[Idx - URes];
//                 NewPos3D[Idx] += 0.25 * Pos3D[Idx + URes];
//             }
//         }

//         // Project back onto the surface
//         #pragma omp parallel for
//         for (int j = 1; j < VRes - 1; ++j)
//         {
//             for (int i = 1; i < URes - 1; ++i)
//             {
//                 int Idx = j * URes + i;
//                 Eigen::Vector2d uv{ U(j, i), V(j, i) };
//                 Eigen::Vector3d p = ImgSampler.Sample(uv, dfy::ImgInterpType::CUBIC);
//                 Eigen::Matrix<double, 3, 2> J;
//                 // Projection by Newton iteration
//                 for (int k = 0; k < 10; ++k)
//                 {
//                     J = ImgSampler.Jacobian(uv, dfy::ImgInterpType::CUBIC);
//                     uv -= (J.transpose() * J).eval().inverse() * (J.transpose() * (p - NewPos3D[Idx])).eval();
//                     p = ImgSampler.Sample(uv, dfy::ImgInterpType::CUBIC);
//                 }
//                 U(j, i) = 0.99 * U(j, i) + 0.01 * uv.x();
//                 V(j, i) = 0.99 * V(j, i) + 0.01 * uv.y();
//                 Pos3D[Idx] = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//             }
//         }

//         // Compute residual
//         double NewErr = 0.0;
//         for (int j = 0; j < VRes; ++j)
//         {
//             for (int i = 0; i < URes; ++i)
//             {
//                 int Idx = j * URes + 1;
//                 Eigen::Vector3d Pij = Pos3D[Idx];
//                 if (i < URes - 1)
//                 {
//                     Eigen::Vector3d Po = Pos3D[Idx + 1];
//                     NewErr += (Pij - Po).squaredNorm();
//                 }
//                 if (j < VRes - 1)
//                 {
//                     Eigen::Vector3d Po = Pos3D[Idx + URes];
//                     NewErr += (Pij - Po).squaredNorm();
//                 }
//             }
//         }
//         if (std::abs(NewErr - Err) < 1e-6)
//         {
//             std::cout << "Converged to " << NewErr << " after " << Iter + 1 << " iterations." << std::endl;
//             break;
//         }
//         Err = NewErr;
//         std::cout << Err << std::endl;
//     }
//     if (Err > 1e-4)
//         std::cerr << "WARNING: Optimization did not converge." << std::endl;
//     std::cout << Err << std::endl;

//     // Create the quad mesh
//     Eigen::MatrixXd Verts(URes * VRes, 3);
//     int CurIdx = 0;
//     for (int j = 0; j < VRes; ++j)
//     {
//         for (int i = 0; i < URes; ++i)
//             Verts.row(CurIdx++) = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//     }

//     if (Weld)
//         throw std::runtime_error("Boundary welding is not yet implemented.");

//     return dfy::QuadMesh(Verts, URes, VRes);
// }

// dfy::QuadMesh dfy::GImage::AsQuadMesh(int URes, int VRes, bool Weld) const
// {
//     // Initialize an image sampler
//     dfy::ImageSampler ImgSampler(*this);

//     // Initialize grid
//     Eigen::MatrixXd U(VRes, URes);
//     Eigen::MatrixXd V(VRes, URes);
//     for (int j = 0; j < VRes; ++j)
//     {
//         double v = j / (VRes - 1.0);
//         for (int i = 0; i < URes; ++i)
//         {
//             double u = i / (URes - 1.0);
//             U(j, i) = u;
//             V(j, i) = v;
//         }
//     }

//     // Map edges to indices
//     std::map<std::pair<int, int>, int> EMap;
//     for (int j = 0; j < VRes; ++j)
//     {
//         for (int i = 0; i < URes; ++i)
//         {
//             int Idx = j * URes + i;
//             std::pair<int, int> e;
//             e.first = Idx;
//             if (i < URes - 1 && j > 0 && j < VRes - 1)
//             {
//                 e.second = Idx + 1;
//                 EMap.emplace(e, (int)EMap.size());
//             }
//             if (j < VRes - 1 && i > 0 && i < URes - 1)
//             {
//                 e.second = Idx + URes;
//                 EMap.emplace(e, (int)EMap.size());
//             }
//             // if (i < URes - 1 && j < VRes - 1)
//             // {
//             //     e.second = Idx + URes + 1;
//             //     EMap.emplace(e, (int)EMap.size());
//             // }
//             // if (i > 0 && j < VRes - 1)
//             // {
//             //     e.second = Idx + URes - 1;
//             //     EMap.emplace(e, (int)EMap.size());
//             // }
//         }
//     }
//     std::cout << EMap.size() << "-by-" << 2 * (URes - 2) * (VRes - 2) << std::endl;
    
//     // Optimize
//     const int MaxIter = 100;
//     double Err = std::numeric_limits<double>::infinity();
//     std::vector<Eigen::Triplet<double>> JTrips;
//     JTrips.reserve(EMap.size() * 4);
//     Eigen::VectorXd Residual;
//     Eigen::SparseMatrix<double> J;
//     J.resize(EMap.size(), 2 * (URes - 2) * (VRes - 2));
//     std::vector<std::pair<int, int>> DeltaEdges = {
//         std::pair<int, int>{ -1,  0 },
//         std::pair<int, int>{  1,  0 },
//         std::pair<int, int>{  0, -1 },
//         std::pair<int, int>{  0,  1 },
//         // std::pair<int, int>{  1,  1 },
//         // std::pair<int, int>{ -1, -1 },
//     };
//     for (int Iter = 0; Iter < MaxIter && Err > 1e-4; ++Iter)
//     {
//         // Create Jacobian and residual vector
//         JTrips.clear();
//         Residual.setZero(EMap.size());

//         // #pragma omp parallel for
//         for (int j = 1; j < VRes - 1; ++j)
//         {
//             for (int i = 1; i < URes - 1; ++i)
//             {
//                 int GlobalIdx = j * URes + i;
//                 int LocalIdx = (j - 1) * (URes - 2) + i - 1;

//                 Eigen::Vector3d Phi = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//                 Eigen::Matrix<double, 3, 2> J = ImgSampler.Jacobian(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//                 for (int k = 0; k < DeltaEdges.size(); ++k)
//                 {
//                     int di, dj;
//                     std::tie(di, dj) = DeltaEdges[k];
//                     Eigen::Vector3d Phi2 = ImgSampler.Sample(U(j + dj, i + di), V(j + dj, i + di), dfy::ImgInterpType::CUBIC);
//                     Eigen::Vector3d Phi12 = Phi - Phi2;
//                     Eigen::Vector2d E12{ U(j, i), V(j, i) };
//                     E12 -= Eigen::Vector2d{ U(j + dj, i + di), V(j + dj, i + di) };
//                     Eigen::Vector2d GE = 2 * (E12 + J.transpose() * Phi12);
//                     // Eigen::Vector2d GE = 2 * E12;
//                     // Eigen::Vector2d GE = 2 * J.transpose() * Phi12;

//                     int OtherIdx = (j + dj) * URes + (i + di);
//                     if (OtherIdx > GlobalIdx)
//                         GE = -GE;
//                     std::pair<int, int> EIdx{ std::min(OtherIdx, GlobalIdx), std::max(OtherIdx, GlobalIdx) };
//                     JTrips.emplace_back(EMap[EIdx], 2 * LocalIdx, GE.x());
//                     JTrips.emplace_back(EMap[EIdx], 2 * LocalIdx + 1, GE.y());
//                     Residual[EMap[EIdx]] = 1e-1 * E12.squaredNorm() + 1e-2 * Phi12.squaredNorm();
//                     // Residual[EMap[EIdx]] = E12.squaredNorm();
//                     // Residual[EMap[EIdx]] = 2.5e-2 * Phi12.squaredNorm();
//                 }
//             }
//         }


//         J.setFromTriplets(JTrips.begin(), JTrips.end());

//         Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> Solver;
//         Solver.compute(J);
//         Eigen::VectorXd UV = Solver.solve(Residual);

//         Err = Residual.squaredNorm();

//         for (int j = 1; j < VRes - 1; ++j)
//         {
//             for (int i = 1; i < URes - 1; ++i)
//             {
//                 int Idx = (j - 1) * (URes - 2) + i - 1;
//                 U(j, i) -= 1.0 / URes * UV[2 * Idx];
//                 V(j, i) -= 1.0 / VRes * UV[2 * Idx + 1];
//             }
//         }
//     }
//     if (Err > 1e-4)
//         std::cerr << "WARNING: Optimization did not converge." << std::endl;
//     std::cout << Err << std::endl;

//     // Create the quad mesh
//     Eigen::MatrixXd Verts(URes * VRes, 3);
//     int CurIdx = 0;
//     for (int j = 0; j < VRes; ++j)
//     {
//         for (int i = 0; i < URes; ++i)
//             Verts.row(CurIdx++) = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//     }

//     if (Weld)
//         throw std::runtime_error("Boundary welding is not yet implemented.");

//     return dfy::QuadMesh(Verts, URes, VRes);
// }

// dfy::QuadMesh dfy::GImage::AsQuadMesh(const dfy::Mesh& M, int URes, int VRes, bool Weld) const
// {
//     // Initialize an image sampler
//     dfy::ImageSampler ImgSampler(*this);

//     // Initialize grid
//     Eigen::MatrixXd U(VRes, URes);
//     Eigen::MatrixXd V(VRes, URes);
//     for (int j = 0; j < VRes; ++j)
//     {
//         double v = j / (VRes - 1.0);
//         for (int i = 0; i < URes; ++i)
//         {
//             double u = i / (URes - 1.0);
//             U(j, i) = u;
//             V(j, i) = v;
//         }
//     }

//     // Iterate until convergence or max iterations
//     const int MaxIter = 10;
//     const int NSamples = (1 << 15) / std::max(URes, VRes);
//     double Err = 0.0;
//     for (int j = 0; j < VRes; ++j)
//     {
//         for (int i = 0; i < URes; ++i)
//         {
//             Eigen::Vector3d Pij = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//             if (i < URes - 1)
//             {
//                 Eigen::Vector3d Po = ImgSampler.Sample(U(j, i + 1), V(j, i + 1), dfy::ImgInterpType::CUBIC);
//                 Err += (Pij - Po).squaredNorm();
//             }
//             if (j < VRes - 1)
//             {
//                 Eigen::Vector3d Po = ImgSampler.Sample(U(j + 1, i), V(j + 1, i), dfy::ImgInterpType::CUBIC);
//                 Err += (Pij - Po).squaredNorm();
//             }
//         }
//     }
//     Eigen::MatrixXd U1 = U;
//     Eigen::MatrixXd U2 = U;
//     Eigen::MatrixXd V1 = V;
//     Eigen::MatrixXd V2 = V;
//     double BestErr = Err;
//     Eigen::MatrixXd BestU = U;
//     Eigen::MatrixXd BestV = V;
//     for (int Iter = 0; Iter < MaxIter; ++Iter)
//     {
//         // Resample each row
//         #pragma omp parallel for
//         for (int j = 0; j < VRes; ++j)
//         {
//             // Allocate memory
//             std::vector<Eigen::Vector2d> Samples;
//             std::vector<double> ArcLength;
//             Samples.reserve(std::max(URes, VRes) * NSamples);
//             ArcLength.reserve(std::max(URes, VRes) * NSamples);
//             // Compute samples
//             Samples.clear();
//             for (int i = 0; i < URes - 1; ++i)
//             {
//                 for (int k = 0; k < NSamples; ++k)
//                 {
//                     double t = k / ((double)NSamples);
//                     double u = (1 - t) * U(j, i) + t * U(j, i + 1);
//                     double v = (1 - t) * V(j, i) + t * V(j, i + 1);
//                     Samples.emplace_back(u, v);
//                 }
//             }
//             Samples.emplace_back(U(j, URes - 1), V(j, URes - 1));

//             // Compute arc length
//             ArcLength.clear();
//             ArcLength.emplace_back(0.0);
//             Eigen::Vector3d LastPos = ImgSampler.Sample(Samples[0], dfy::ImgInterpType::CUBIC);
//             Eigen::Vector3d CurPos = LastPos;
//             for (int k = 1; k < Samples.size(); ++k)
//             {
//                 CurPos = ImgSampler.Sample(Samples[k], dfy::ImgInterpType::CUBIC);
//                 ArcLength.emplace_back(ArcLength[k - 1] + (CurPos - LastPos).eval().norm());
//                 LastPos = CurPos;
//             }

//             // Resample according to uniform arc length
//             double SampleStep = ArcLength.back() / URes;
//             double NextLength = SampleStep;
//             int CurCol = 1;
//             for (int k = 1; k < ArcLength.size(); ++k)
//             {
//                 // First and last positions are kept fixed
//                 if (CurCol == (URes - 1))
//                     break;

//                 // Skip samples until we reach a suitable arc length
//                 if (ArcLength[k] < NextLength)
//                     continue;
                
//                 // Use this position as the current sample
//                 U1(j, CurCol) = Samples[k].x();
//                 V1(j, CurCol) = Samples[k].y();
//                 // Update the length for the next sample
//                 NextLength += SampleStep;
//                 // Update the current column
//                 CurCol++;
//             }
//         }


//         // Resample each column
//         #pragma omp parallel for
//         for (int i = 0; i < URes; ++i)
//         {
//             // Allocate memory
//             std::vector<Eigen::Vector2d> Samples;
//             std::vector<double> ArcLength;
//             Samples.reserve(std::max(URes, VRes) * NSamples);
//             ArcLength.reserve(std::max(URes, VRes) * NSamples);
//             // Compute samples
//             Samples.clear();
//             for (int j = 0; j < VRes - 1; ++j)
//             {
//                 for (int k = 0; k < NSamples; ++k)
//                 {
//                     double t = k / ((double)NSamples);
//                     double u = (1 - t) * U(j, i) + t * U(j + 1, i);
//                     double v = (1 - t) * V(j, i) + t * V(j + 1, i);
//                     Samples.emplace_back(u, v);
//                 }
//             }
//             Samples.emplace_back(U(VRes - 1, i), V(VRes - 1, i));

//             // Compute arc length
//             ArcLength.clear();
//             ArcLength.emplace_back(0.0);
//             Eigen::Vector3d LastPos = ImgSampler.Sample(Samples[0], dfy::ImgInterpType::CUBIC);
//             Eigen::Vector3d CurPos = LastPos;
//             for (int k = 1; k < Samples.size(); ++k)
//             {
//                 CurPos = ImgSampler.Sample(Samples[k], dfy::ImgInterpType::CUBIC);
//                 ArcLength.emplace_back(ArcLength[k - 1] + (CurPos - LastPos).eval().norm());
//                 LastPos = CurPos;
//             }

//             // Resample according to uniform arc length
//             double SampleStep = ArcLength.back() / VRes;
//             double NextLength = SampleStep;
//             int CurRow = 1;
//             for (int k = 1; k < ArcLength.size(); ++k)
//             {
//                 // First and last positions are kept fixed
//                 if (CurRow == (VRes - 1))
//                     break;

//                 // Skip samples until we reach a suitable arc length
//                 if (ArcLength[k] < NextLength)
//                     continue;
                
//                 // Use this position as the current sample
//                 U2(CurRow, i) = Samples[k].x();
//                 V2(CurRow, i) = Samples[k].y();
//                 // Update the length for the next sample
//                 NextLength += SampleStep;
//                 // Update the current row
//                 CurRow++;
//             }
//         }

//         // Average rows optimization and column optimization
//         U = 0.5 * (U1 + U2);
//         V = 0.5 * (V1 + V2);

//         // Compute residual
//         double NewErr = 0.0;
//         for (int j = 0; j < VRes; ++j)
//         {
//             for (int i = 0; i < URes; ++i)
//             {
//                 Eigen::Vector3d Pij = ImgSampler.Sample(U(j, i), V(j, i), dfy::ImgInterpType::CUBIC);
//                 if (i < URes - 1)
//                 {
//                     Eigen::Vector3d Po = ImgSampler.Sample(U(j, i + 1), V(j, i + 1), dfy::ImgInterpType::CUBIC);
//                     NewErr += (Pij - Po).squaredNorm();
//                 }
//                 if (j < VRes - 1)
//                 {
//                     Eigen::Vector3d Po = ImgSampler.Sample(U(j + 1, i), V(j + 1, i), dfy::ImgInterpType::CUBIC);
//                     NewErr += (Pij - Po).squaredNorm();
//                 }
//             }
//         }
//         if (std::abs(NewErr - Err) < 1e-6)
//         {
//             std::cout << "Converged to " << NewErr << " after " << Iter + 1 << " iterations." << std::endl;
//             break;
//         }
//         Err = NewErr;
//         if (Err < BestErr)
//         {
//             BestU = U;
//             BestV = V;
//             BestErr = Err;
//         }
//         std::cout << Err << " (best " << BestErr << ")" << std::endl;
//     }

//     // Create the quad mesh
//     Eigen::MatrixXd Verts(URes * VRes, 3);
//     int CurIdx = 0;
//     for (int j = 0; j < VRes; ++j)
//     {
//         for (int i = 0; i < URes; ++i)
//             Verts.row(CurIdx++) = ImgSampler.Sample(BestU(j, i), BestV(j, i), dfy::ImgInterpType::CUBIC);
//             // Verts.row(CurIdx++) = Eigen::Vector3d{ BestU(j, i), BestV(j, i), 0};
//     }

//     if (Weld)
//         throw std::runtime_error("Boundary welding is not yet implemented.");

//     return dfy::QuadMesh(Verts, URes, VRes);
// }