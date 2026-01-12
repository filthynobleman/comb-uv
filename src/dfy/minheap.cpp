/**
 * @file        minheap.cpp
 * 
 * @brief       Implements dfy::MinHeap.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-25
 */
#include <dfy/minheap.hpp>
#include <cmath>
#include <iostream>

std::pair<double, size_t> dfy::MinHeap::FindMin() const
{
    return { m_Sign * m_Nodes[0].first, m_Nodes[0].second };
}

size_t dfy::MinHeap::Size() const { return m_Nodes.size(); }


void dfy::MinHeap::MoveUp(size_t Element)
{
    size_t v = m_Perm[Element];
    while (v != 0)
    {
        size_t p = (v - 1) >> 1;  // parent is floor(v / 2)
        // if parent is already smaller, it's done
        if (m_Nodes[p].first < m_Nodes[v].first)
            break;

        // exchange node with parent
        std::swap(m_Nodes[p], m_Nodes[v]);

        // remember to swap the permutation to preserve the pointers
        std::swap(m_Perm[m_Nodes[p].second], m_Perm[m_Nodes[v].second]);

        // Next iteration
        v = p;
    }
}

void dfy::MinHeap::MoveDown(size_t Element)
{
    size_t v = m_Perm[Element];
    while (true)
    {
        size_t u = (v << 1) + 1; // left child is 2 * v
        // If no children, break
        if (u >= Size())
            break;
        // If two children, pick minimum key
        if (u < Size() - 1)
        {
            if (m_Nodes[u].first > m_Nodes[u + 1].first)
                u++;
        }

        // If children has a larger key, break
        if (m_Nodes[u].first > m_Nodes[v].first)
            break;
        
        // Swap node and children
        std::swap(m_Nodes[v], m_Nodes[u]);

        // remember to swap the permutation to preserve the pointers
        std::swap(m_Perm[m_Nodes[v].second], m_Perm[m_Nodes[u].second]);

        // Next iteration
        v = u;
    }
}

double dfy::MinHeap::GetKey(size_t Element) const
{
    return m_Sign * m_Nodes[m_Perm[Element]].first;
}

void dfy::MinHeap::DecreaseKey(size_t Element, double Decrement)
{
    size_t v = m_Perm[Element];
    m_Nodes[v].first -= m_Sign * Decrement;

    // If min-heap, keys are stored positive, decrement reduces the key
    // value and stuff must move up
    if (m_Sign > 0)
        MoveUp(Element);
    // If max-heap, keys are stored negative (but shown positive), decrement
    // increases the key value, and stuff moves down
    else
        MoveDown(Element);
}

void dfy::MinHeap::IncreaseKey(size_t Element, double Increment)
{
    size_t v = m_Perm[Element];
    m_Nodes[v].first += m_Sign * Increment;

    // If min-heap, keys are stored positive, increment increases the
    // key value, and stuf moves down
    if (m_Sign > 0)
        MoveDown(Element);
    // If max-heap, keys are stored negative (but shown positive), increment
    // decreases the key value, and stuff moves up
    else
        MoveUp(Element);
}

void dfy::MinHeap::SetKey(size_t Element, double NewKey)
{
    NewKey *= m_Sign;
    size_t v = m_Perm[Element];
    // If the operation increases the key, move down
    if (NewKey > m_Nodes[v].first)
    {
        m_Nodes[v].first = NewKey;
        MoveDown(Element);
    }
    // If the operation decreases the key, move up
    else if (NewKey < m_Nodes[v].first)
    {
        m_Nodes[v].first = NewKey;
        MoveUp(Element);
    }
    // If key is unchanged, there is no need to modify the heap
}


void dfy::MinHeap::Insert(double Key)
{
    m_Perm.push_back(m_Nodes.size());
    m_Nodes.emplace_back(m_Sign * Key, m_Nodes.size());
    MoveUp(m_Nodes.size() - 1);
}


dfy::MinHeap::MinHeap(const std::vector<double>& Keys,
                      bool IsMaxHeap)
{
    // Min-heap store keys with the right sign.
    // Max-heap store keys with negated sign.
    m_Sign = 1;
    if (IsMaxHeap)
        m_Sign = -1;
    
    // Add one by one and update the heap
    for (double k : Keys)
        Insert(k);
}

dfy::MinHeap::MinHeap(const double* const Keys,
                      size_t NumKeys,
                      bool IsMaxHeap)
{
    // Min-heap store keys with the right sign.
    // Max-heap store keys with negated sign.
    m_Sign = 1;
    if (IsMaxHeap)
        m_Sign = -1;
    
    // Add one by one and update the heap
    for (size_t i = 0; i < NumKeys; ++i)
        Insert(Keys[i]);
}

dfy::MinHeap::~MinHeap() { }