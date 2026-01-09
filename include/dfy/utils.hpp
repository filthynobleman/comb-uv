/**
 * @file        utils.hpp
 * 
 * @brief       
 * 
 * @author      Filippo Maggioli (maggioli.filippo@gmail.com)
 * 
 * @date        2026-01-09
 */
#pragma once

#include <utility>
#include <functional>

#include <Eigen/Eigen>

#define NOMINMAX


namespace rmt
{

template<typename T>
struct PairHash
{
    std::size_t operator()(const std::pair<T, T>& X) const noexcept
    {
        std::hash<T> h{};
        return h(X.first) ^ h(X.second);
    }
};

template<typename T>
struct TripleHash
{
    std::size_t operator()(const std::tuple<T, T, T>& X) const noexcept
    {
        std::hash<T> h{};
        return h(std::get<0>(X)) ^ (h(std::get<1>(X) ^ (h(std::get<2>(X)) << 1)) << 1);
    }
};




} // namespace rmt

namespace Eigen
{
    
inline Eigen::internal::all_t all = Eigen::placeholders::all;

} // namespace Eigen