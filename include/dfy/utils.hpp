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
#include <cassert>

#include <Eigen/Eigen>

#define NOMINMAX

namespace dfy
{
    
enum UVMapAlgorithm
{
    TUTTE,
    HARMONIC,
    ARAP,
    CONFORMAL
};

} // namespace dfy


namespace Eigen
{
    
inline Eigen::internal::all_t all = Eigen::placeholders::all;

} // namespace Eigen