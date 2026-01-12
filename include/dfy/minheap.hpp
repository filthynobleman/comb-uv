/**
 * @file        minheap.hpp
 * 
 * @brief       A class for implementing a min-heap data structure.
 * 
 * @details     This file contains the declaration of the class dfy::MinHeap,
 *              which provide an interface for a min-heap data structure.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-25
 */
#pragma once

#include <vector>
#include <cstddef>

namespace dfy
{

/**
 * @brief       This class represents a min-heap data structure.
 * 
 * @details     The class dfy::MinHeap provides an interface for a min-heap
 *              data structure.\n 
 *              The dfy::MinHeap contains integer elements mapped to double
 *              precision floating point keys.\n 
 *              While this structure seems very specific, it can be used for
 *              basically everything, since key updates can be deferred outside
 *              of this class, and integer values can be used to map array
 *              positions.\n 
 *              This class can also be used to realize max-heap structures
 *              (see constructor's documentation). 
 * 
 * @warning     This interface does not support deletions or insertions,
 *              but only increasing and decreasing of keys.\n 
 *              This limitation is intentional, as it allows for an efficient
 *              implementation.
 */
class MinHeap
{
private:
    /**
     * @brief       Heap elements.
     * 
     * @details     This vector contains the heap elements. Specifically, each
     *              element in this vector determines the key associated to a
     *              given element, and the vector itself represents the heap tree.
     */
    std::vector<std::pair<double, size_t>> m_Nodes;

    /**
     * @brief       The inverse permutation of the tree.
     * 
     * @details     This vector contains the inverse permutation of the heap tree.\n 
     *              Namely, the value \code m_Perms[i] \endcode tells the position
     *              in the heap-tree where the element \code i \endcode is stored.
     */
    std::vector<size_t> m_Perm;

    /**
     * @brief       The comparison sign.
     * 
     * @details     This is the sign for handling the comparison.\n 
     *              If this value is positive, the dfy::MinHeap implements
     *              a min-heap. If this value is negative, the dfy::MinHeap
     *              implements a max-heap.
     */
    int m_Sign;


    void MoveUp(size_t Element);
    void MoveDown(size_t Element);
    void Insert(double Key);


public:
    /**
     * @brief       Create a new dfy::MinHeap from the vector of keys.
     * 
     * @details     This constructor creates a new dfy::MinHeap from the
     *              given vector of keys.\n 
     *              The key \code Keys[i] \endcode from the provided vector of keys
     *              is associated to the integer value \code i \endcode before
     *              the heapifying of the structure.\n 
     *              If the optional parameter \code AsMaxHeap \endcode is set to
     *              true (default is false), the dfy::MinHeap is modified to
     *              implement a max-heap data structure.
     * 
     * @param Keys The vector of keys to initialize the heap.
     * @param AsMaxHeap Decide whether or not the heap must implement a min- or max- heap.
     */
    MinHeap(const std::vector<double>& Keys,
            bool AsMaxHeap = false);

    /**
     * @brief       Create a new dfy::MinHeap from the array of keys.
     * 
     * @details     This constructor creates a new dfy::MinHeap from the
     *              given array of keys.\n 
     *              The key \code Keys[i] \endcode from the provided array of keys
     *              is associated to the integer value \code i \endcode before
     *              the heapifying of the structure.\n 
     *              If the optional parameter \code AsMaxHeap \endcode is set to
     *              true (default is false), the dfy::MinHeap is modified to
     *              implement a max-heap data structure.
     * 
     * @param Keys The array of keys to initialize the heap.
     * @param NumKeys The number of keys in the array.
     * @param AsMaxHeap Decide whether or not the heap must implement a min- or max- heap.
     */
    MinHeap(const double* const Keys,
            size_t NumKeys,
            bool AsMaxHeap = false);

    /**
     * @brief       Default destructor.
     * 
     * @details     Default destructor.
     */
    ~MinHeap();


    /**
     * @brief       Returns the size of the min-heap.
     * 
     * @details     This method returns the size of the min-heap, that is
     *              the number of elements contained in the heap.
     * 
     * @return size_t The heap's size.
     */
    size_t Size() const;


    /**
     * @brief       Returns the minimum (or maximum) element.
     * 
     * @details     This method returns the element of the heap with the
     *              minimum key.\n 
     *              If the heap is a max-heap, this method returns the element
     *              with maximum key.
     * 
     * @return std::pair<double, size_t> The kay-value pair in the heap with minimum key.
     * 
     * @throws dfy::OutOfBoundError if the heap is empty.
     */
    std::pair<double, size_t> FindMin() const;

    /**
     * @brief       Returns the key of the element.
     * 
     * @details     This method returns the key associated to the given
     *              element.
     * 
     * @param Element An element of the heap.
     * @return double The key of the element.
     * 
     * @throws dfy::OutOfBoundError if Element is larger than the heap size.
     */
    double GetKey(size_t Element) const;

    /**
     * @brief       Decreases the key of the element by the given value.
     * 
     * @details     This method decreases the key associated with the given
     *              element by the given value, and then restores the heap.\n 
     *              This operation is a decrement independently if the heap
     *              is a min-heap or a max-heap.
     * 
     * @param Element The element whose key must be decremented.
     * @param Decrement The decrement to operate on the key.
     * 
     * @throws dfy::OutOfBoundError if Element is larger than the heap size.
     */
    void DecreaseKey(size_t Element,
                     double Decrement);

    /**
     * @brief       Increases the key of the element by the given value.
     * 
     * @details     This method increases the key associated with the given
     *              element by the given value, and then restores the heap.\n 
     *              This operation is an increment independently if the heap
     *              is a min-heap or a max-heap.
     * 
     * @param Element The element whose key must be incremented.
     * @param Increment The increment to operate on the key.
     * 
     * @throws dfy::OutOfBoundError if Element is larger than the heap size.
     */
    void IncreaseKey(size_t Element,
                     double Increment);

    /**
     * @brief       Set the key of the element to a given value.
     * 
     * @details     This method updates the key associated to the given element
     *              with the given value, and then restores the heap.
     * 
     * @param Element The element whose key must be modified.
     * @param NewKey The new value of the key.
     */
    void SetKey(size_t Element,
                double NewKey);
};

} // namespace dfy