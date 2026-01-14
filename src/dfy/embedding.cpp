/**
 * @file        embedding.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-13
 */
#include <dfy/embedding.hpp>

#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>

#include <dfy/io.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

dfy::Embedding::Embedding(const dfy::Mesh &M)
    : m_Mesh(M)
{
    igl::boundary_loop(GetMesh().Triangles(), m_BLoop);
    m_UV.setZero(m_Mesh.NumVertices(), 2);
}

dfy::Embedding::Embedding(const dfy::Embedding &E)
    : m_Mesh(E.m_Mesh)
{
    m_BLoop = E.m_BLoop;
    m_UV = E.m_UV;
}

dfy::Embedding::Embedding(dfy::Embedding &&E)
    : m_Mesh(E.m_Mesh)
{
    m_BLoop = std::move(E.m_BLoop);
    m_UV = std::move(E.m_UV);
}

dfy::Embedding::~Embedding() { }

const dfy::Mesh &dfy::Embedding::GetMesh() const { return m_Mesh; }
std::vector<int> &dfy::Embedding::GetBoundaryLoop() { return m_BLoop; }
const std::vector<int> &dfy::Embedding::GetBoundaryLoop() const { return m_BLoop; }
Eigen::MatrixXd &dfy::Embedding::GetUV() { return m_UV; }
const Eigen::MatrixXd &dfy::Embedding::UV() const { return m_UV; }

void Circle2Square(Eigen::MatrixXd& P)
{
    int n = P.rows();
    Eigen::VectorXd Theta;
    Theta.resize(n);
    Eigen::Vector4d CornerErr;
    CornerErr.setConstant(std::numeric_limits<double>::infinity());
    Eigen::Vector4i CornerIdx;
    CornerIdx.setConstant(-1);
    for (int i = 0; i < n; ++i)
    {
        Theta[i] = std::atan2(P(i, 1), P(i, 0)) + M_PI;
        for (int j = 0; j < 4; ++j)
        {
            double ThetaStar = M_PI_4 + j * M_PI_2;
            double Err = std::abs(Theta[i] - ThetaStar);
            if (Err < CornerErr[j])
            {
                CornerErr[j] = Err;
                CornerIdx[j] = i;
            }
        }
    }
    P.row(CornerIdx[0]) = Eigen::RowVector2d{  1.0,  1.0 };
    P.row(CornerIdx[1]) = Eigen::RowVector2d{ -1.0,  1.0 };
    P.row(CornerIdx[2]) = Eigen::RowVector2d{ -1.0, -1.0 };
    P.row(CornerIdx[3]) = Eigen::RowVector2d{  1.0, -1.0 };

    for (int j0 = 0; j0 < 4; ++j0)
    {
        int j1 = (j0 + 1) % 4;
        int c0 = CornerIdx[j0];
        int c1 = CornerIdx[j1];
        if (c0 == c1)
            continue;
        for (int i = (c0 + 1) % n; i != c1; i = (i + 1) % n)
        {
            double t = (Theta[i] - Theta[c0]) / (Theta[c1] - Theta[c0]);
            if (j0 == 3 && Theta[i] < M_PI)
                t = (Theta[i] - (Theta[c0] - 2 * M_PI)) / (Theta[c1] - (Theta[c0] - 2 * M_PI));
            else if (j0 == 3)
                t = (Theta[i] - Theta[c0]) / (Theta[c1] + 2 * M_PI - Theta[c0]);
            
            P.row(i) = (1 - t) * P.row(c0) + t * P.row(c1);
        }
    }
}

void Circle2Tri(Eigen::MatrixXd& P)
{
    int n = P.rows();
    Eigen::VectorXd Theta;
    Theta.resize(n);
    Eigen::Vector3d CornerErr;
    CornerErr.setConstant(std::numeric_limits<double>::infinity());
    Eigen::Vector3i CornerIdx;
    CornerIdx.setConstant(-1);
    for (int i = 0; i < n; ++i)
    {
        Theta[i] = std::atan2(P(i, 1), P(i, 0)) + M_PI;
        for (int j = 0; j < 3; ++j)
        {
            double ThetaStar = M_PI_2 + 2.0 * j * M_PI / 3.0;
            double Err = std::abs(Theta[i] - ThetaStar);
            if (Err < CornerErr[j])
            {
                CornerErr[j] = Err;
                CornerIdx[j] = i;
            }
        }
    }
    P.row(CornerIdx[0]) = Eigen::RowVector2d{  0.0,  1.0 };
    P.row(CornerIdx[1]) = Eigen::RowVector2d{ -1.0, -1.0 };
    P.row(CornerIdx[2]) = Eigen::RowVector2d{  1.0, -1.0 };

    for (int j0 = 0; j0 < 3; ++j0)
    {
        int j1 = (j0 + 1) % 3;
        int c0 = CornerIdx[j0];
        int c1 = CornerIdx[j1];
        if (c0 == c1)
            continue;
        for (int i = (c0 + 1) % n; i != c1; i = (i + 1) % n)
        {
            double t = (Theta[i] - Theta[c0]) / (Theta[c1] - Theta[c0]);
            if (j0 == 2 && Theta[i] < M_PI)
                t = (Theta[i] - (Theta[c0] - 2 * M_PI)) / (Theta[c1] - (Theta[c0] - 2 * M_PI));
            else if (j0 == 2)
                t = (Theta[i] - Theta[c0]) / (Theta[c1] + 2 * M_PI - Theta[c0]);
            
            P.row(i) = (1 - t) * P.row(c0) + t * P.row(c1);
        }
    }
}

void Circle2Hex(Eigen::MatrixXd& P)
{
    int n = P.rows();
    Eigen::VectorXd Theta;
    Theta.resize(n);
    Eigen::Vector<double, 6> CornerErr;
    CornerErr.setConstant(std::numeric_limits<double>::infinity());
    Eigen::Vector<int, 6> CornerIdx;
    CornerIdx.setConstant(-1);
    for (int i = 0; i < n; ++i)
    {
        Theta[i] = std::atan2(P(i, 1), P(i, 0)) + M_PI;
        for (int j = 0; j < 6; ++j)
        {
            double ThetaStar = j * M_PI / 3.0;
            double Err = std::abs(Theta[i] - ThetaStar);
            if (Err < CornerErr[j])
            {
                CornerErr[j] = Err;
                CornerIdx[j] = i;
            }
        }
    }
    P.row(CornerIdx[0]) = Eigen::RowVector2d{  1.0,  0.0 };
    P.row(CornerIdx[1]) = Eigen::RowVector2d{  0.5,  1.0 };
    P.row(CornerIdx[2]) = Eigen::RowVector2d{ -0.5,  1.0 };
    P.row(CornerIdx[3]) = Eigen::RowVector2d{ -1.0,  0.0 };
    P.row(CornerIdx[4]) = Eigen::RowVector2d{ -0.5, -1.0 };
    P.row(CornerIdx[5]) = Eigen::RowVector2d{  0.5, -1.0 };

    for (int j0 = 0; j0 < 6; ++j0)
    {
        int j1 = (j0 + 1) % 6;
        int c0 = CornerIdx[j0];
        int c1 = CornerIdx[j1];
        if (c0 == c1)
            continue;
        for (int i = (c0 + 1) % n; i != c1; i = (i + 1) % n)
        {
            double t = (Theta[i] - Theta[c0]) / (Theta[c1] - Theta[c0]);
            if (j0 == 5)
                t = (Theta[i] - Theta[c0]) / (Theta[c1] + 2 * M_PI - Theta[c0]);
            
            P.row(i) = (1 - t) * P.row(c0) + t * P.row(c1);
        }
    }
}

void dfy::Embedding::MapBoundary(dfy::BoundaryMap BMap)
{
    Eigen::MatrixXd BUV;
    Eigen::VectorXi BLoop = Eigen::VectorXi::Map(GetBoundaryLoop().data(), GetBoundaryLoop().size());
    igl::map_vertices_to_circle(GetMesh().Vertices(), BLoop, BUV);
    switch (BMap)
    {
        case dfy::BoundaryMap::SQUARE:
        Circle2Square(BUV);
        break;
        
        case dfy::BoundaryMap::TRIANGLE:
        Circle2Tri(BUV);
        break;
        
        case dfy::BoundaryMap::HEXAGON:
        Circle2Hex(BUV);
        break;
        
        case dfy::BoundaryMap::CIRCLE:
        break;
        
        default:
        throw std::runtime_error("Unknown boundary condition.");
        break;
    }
    dfy::ExportPointCloud("./boundary.obj", BUV);

    GetUV()(BLoop, Eigen::all) = 0.5 * (BUV.array() + 1);
}