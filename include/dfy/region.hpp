/**
 * @file        region.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-12
 */
#pragma once

#include <dfy/manifoldmesh.hpp>

#include <set>

namespace dfy
{
    
class Region
{
private:
    std::set<int> m_Verts;
    std::set<int> m_Edges;
    std::set<int> m_Faces;

public:
    Region();
    Region(const dfy::Region& R);
    Region(dfy::Region&& R);
    dfy::Region& operator=(const dfy::Region& R);
    dfy::Region& operator=(dfy::Region&& R);
    ~Region();

    void Vertices(std::vector<int>& VertsVec) const;
    void Edges(std::vector<int>& EdgesVec) const;
    void Faces(std::vector<int>& FacesVec) const;

    int NumVertices() const;
    int NumEdges() const;
    int NumFaces() const;

    int EulerCharacteristic() const;

    bool IsDisk() const;

    void AddVertex(int VertID);
    void AddEdge(int EdgeID);
    void AddFace(int FaceID);

    friend dfy::Region Union(const dfy::Region&, const dfy::Region);
    friend dfy::Region Intersection(const dfy::Region&, const dfy::Region);
    friend dfy::Region LineIntersection(const dfy::Region&, const dfy::Region);
};

dfy::Region Union(const dfy::Region& R1, const dfy::Region R2);
dfy::Region Intersection(const dfy::Region& R1, const dfy::Region R2);
dfy::Region LineIntersection(const dfy::Region& R1, const dfy::Region R2);


} // namespace dfy
