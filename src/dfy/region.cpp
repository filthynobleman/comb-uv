/**
 * @file        region.cpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#include <dfy/region.hpp>
#include <dfy/utils.hpp>


dfy::Region::Region() 
{
    m_Verts.clear();
    m_Edges.clear();
    m_Faces.clear();
}

dfy::Region::Region(const dfy::Region& R)
{
    m_Verts = R.m_Verts;
    m_Edges = R.m_Edges;
    m_Faces = R.m_Faces;
}

dfy::Region::Region(dfy::Region&& R)
{
    m_Verts = std::move(R.m_Verts);
    m_Edges = std::move(R.m_Edges);
    m_Faces = std::move(R.m_Faces);
}

dfy::Region& dfy::Region::operator=(const dfy::Region& R)
{
    m_Verts = R.m_Verts;
    m_Edges = R.m_Edges;
    m_Faces = R.m_Faces;
    return *this;
}

dfy::Region& dfy::Region::operator=(dfy::Region&& R)
{
    m_Verts = std::move(R.m_Verts);
    m_Edges = std::move(R.m_Edges);
    m_Faces = std::move(R.m_Faces);
    return *this;
}

dfy::Region::~Region() { }

void dfy::Region::Vertices(std::vector<int> &VertsVec) const
{
    VertsVec.clear();
    if (NumVertices() == 0)
        return;
    VertsVec.reserve(m_Verts.size());
    VertsVec.insert(VertsVec.begin(), m_Verts.begin(), m_Verts.end());
    std::sort(VertsVec.begin(), VertsVec.end());
}

void dfy::Region::Edges(std::vector<int> &EdgesVec) const
{
    EdgesVec.clear();
    if (NumEdges() == 0)
        return;
    EdgesVec.reserve(m_Edges.size());
    EdgesVec.insert(EdgesVec.begin(), m_Edges.begin(), m_Edges.end());
    std::sort(EdgesVec.begin(), EdgesVec.end());
}

void dfy::Region::Faces(std::vector<int> &FacesVec) const
{
    FacesVec.clear();
    if (NumFaces() == 0)
        return;
    FacesVec.reserve(m_Faces.size());
    FacesVec.insert(FacesVec.begin(), m_Faces.begin(), m_Faces.end());
    std::sort(FacesVec.begin(), FacesVec.end());
}

int dfy::Region::NumVertices() const { return m_Verts.size(); }
int dfy::Region::NumEdges() const { return m_Edges.size(); }
int dfy::Region::NumFaces() const { return m_Faces.size(); }

int dfy::Region::EulerCharacteristic() const
{
    return NumVertices() + NumFaces() - NumEdges();
}

bool dfy::Region::IsDisk() const
{
    return EulerCharacteristic() == 1;
}

bool dfy::Region::IsEmpty() const
{
    return (NumVertices() + NumEdges() + NumFaces()) == 0;
}

void dfy::Region::AddVertex(int VertID) { m_Verts.insert(VertID); }
void dfy::Region::AddEdge(int EdgeID) { m_Edges.insert(EdgeID); }
void dfy::Region::AddFace(int FaceID) { m_Faces.insert(FaceID); }

void dfy::Region::Reserve(int nV, int nE, int nT)
{
    m_Verts.reserve(nV);
    m_Edges.reserve(nE);
    m_Faces.reserve(nT);
}
void dfy::Region::Clear()
{
    m_Verts.clear();
    m_Edges.clear();
    m_Faces.clear();
}

dfy::Region dfy::Union(const dfy::Region &R1, const dfy::Region& R2)
{
    dfy::Region R(R1);
    R.Reserve(R.NumVertices() + R2.NumVertices(),
              R.NumEdges() + R2.NumEdges(),
              R.NumFaces() + R2.NumFaces());
    for (int v : R2.m_Verts)
        R.AddVertex(v);
    for (int e : R2.m_Edges)
        R.AddEdge(e);
    for (int f : R2.m_Faces)
        R.AddFace(f);
        
    // std::set_union(R1.m_Verts.begin(), R1.m_Verts.end(),
    //                R2.m_Verts.begin(), R2.m_Verts.end(),
    //                std::inserter(R.m_Verts, R.m_Verts.begin()));
    // std::set_union(R1.m_Edges.begin(), R1.m_Edges.end(),
    //                R2.m_Edges.begin(), R2.m_Edges.end(),
    //                std::inserter(R.m_Edges, R.m_Edges.begin()));
    // std::set_union(R1.m_Faces.begin(), R1.m_Faces.end(),
    //                R2.m_Faces.begin(), R2.m_Faces.end(),
    //                std::inserter(R.m_Faces, R.m_Faces.begin()));
    return R;
}

dfy::Region dfy::Intersection(const dfy::Region &R1, const dfy::Region& R2)
{
    dfy::Region R;
    R.Reserve(R1.NumVertices(), R1.NumEdges(), R1.NumFaces());
    for (int v : R1.m_Verts)
    {
        if (R2.m_Verts.find(v) != R2.m_Verts.end())
            R.AddVertex(v);
    }
    for (int e : R1.m_Edges)
    {
        if (R2.m_Edges.find(e) != R2.m_Edges.end())
            R.AddEdge(e);
    }
    for (int f : R1.m_Faces)
    {
        if (R2.m_Faces.find(f) != R2.m_Faces.end())
            R.AddFace(f);
    }
    // std::set_intersection(R1.m_Verts.begin(), R1.m_Verts.end(),
    //                       R2.m_Verts.begin(), R2.m_Verts.end(),
    //                       std::inserter(R.m_Verts, R.m_Verts.begin()));
    // std::set_intersection(R1.m_Edges.begin(), R1.m_Edges.end(),
    //                       R2.m_Edges.begin(), R2.m_Edges.end(),
    //                       std::inserter(R.m_Edges, R.m_Edges.begin()));
    // std::set_intersection(R1.m_Faces.begin(), R1.m_Faces.end(),
    //                       R2.m_Faces.begin(), R2.m_Faces.end(),
    //                       std::inserter(R.m_Faces, R.m_Faces.begin()));
    return R;
}

dfy::Region dfy::LineIntersection(const dfy::Region &R1, const dfy::Region& R2)
{
    dfy::Region R;
    for (int v : R1.m_Verts)
    {
        if (R2.m_Verts.find(v) != R2.m_Verts.end())
            R.AddVertex(v);
    }
    for (int e : R1.m_Edges)
    {
        if (R2.m_Edges.find(e) != R2.m_Edges.end())
            R.AddEdge(e);
    }
    // std::set_intersection(R1.m_Verts.begin(), R1.m_Verts.end(),
    //                       R2.m_Verts.begin(), R2.m_Verts.end(),
    //                       std::inserter(R.m_Verts, R.m_Verts.begin()));
    // std::set_intersection(R1.m_Edges.begin(), R1.m_Edges.end(),
    //                       R2.m_Edges.begin(), R2.m_Edges.end(),
    //                       std::inserter(R.m_Edges, R.m_Edges.begin()));
    return R;
}