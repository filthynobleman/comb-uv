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
#include <igl/writeOBJ.h>

#include <stb_image_write.h>


bool dfy::ExportMesh(const std::string &Filename, 
                     const dfy::Mesh &M)
{
    return igl::write_triangle_mesh(Filename, M.Vertices(), M.Triangles());
}

bool dfy::ExportMesh(const std::string &Filename, 
                     const dfy::Mesh &M, 
                     const Eigen::MatrixXd &UV, 
                     const Eigen::MatrixXi &TUV,
                     bool Smooth)
{
    if (Smooth)
        return igl::writeOBJ(Filename, 
                             M.Vertices(), M.Triangles(), 
                             M.VertNormals(), M.Triangles(), 
                             UV, TUV);
    else
    {
        int NTris = M.NumTriangles();
        Eigen::MatrixXi FN = Eigen::VectorXi::LinSpaced(NTris, 0, NTris - 1).replicate(1, 3);
        return igl::writeOBJ(Filename, 
                             M.Vertices(), M.Triangles(), 
                             M.FaceNormals(), FN, 
                             UV, TUV);
    }
}

bool dfy::ExportPointCloud(const std::string &Filename, 
                           const Eigen::MatrixXd &Points)
{
    std::ofstream Stream;
    Stream.open(Filename, std::ios::out);
    if (!Stream.is_open())
        return false;

    Stream << "o " << Filename.substr(0, Filename.rfind('.')) << '\n';

    for (int i = 0; i < Points.rows(); ++i)
    {
        Stream << "v ";
        for (int j = 0; j < Points.cols(); ++j)
            Stream << Points(i, j) << ' ';
        for (int j = Points.cols(); j < 3; ++j)
            Stream << 0 << ' ';
        Stream << '\n';
    }

    Stream.close();
    return true;
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

bool dfy::ExportGImage(const std::string &Filename, 
                       const dfy::GImage &Img)
{
    size_t Stride = Img.GetWidth() * 3 * sizeof(unsigned char);
    size_t BufSize = Stride * Img.GetHeight();
    unsigned char* imgdata = (unsigned char*)std::malloc(BufSize);
    if (imgdata == nullptr)
        return false;
    
    for (int j = 0; j < Img.GetHeight(); ++j)
    {
        for (int i = 0; i < Img.GetWidth(); ++i)
        {
            size_t PIdx = (j * Img.GetWidth() + i) * 3;
            Img.GetPixel(i, j, imgdata[PIdx], imgdata[PIdx + 1], imgdata[PIdx + 2]);
        }
    }

    stbi_flip_vertically_on_write(true);
    int Status = stbi_write_png(Filename.c_str(),
                                Img.GetWidth(), Img.GetHeight(), 3,
                                imgdata, Stride);

    delete imgdata;
    return Status != 0;
}