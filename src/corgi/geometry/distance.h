#pragma once

#include <tuple>
#include <cmath>

namespace corgi { 
  
  namespace geom {

    // L_1 norm
    template<std::size_t D, typename T>
    auto manhattan_distance(const std::array<T, D> arr) {
      return [&]<std::size_t... I>(std::index_sequence<I...>) {
        return (std::abs(arr[I]) + ...);
      }(std::make_index_sequence<D>());
    }

    // L_2 norm
    template<std::size_t D, typename T>
    auto eulerian_distance(const std::array<T, D> arr) {
      return [&]<std::size_t... I>(std::index_sequence<I...>) {
        return std::sqrt((std::pow(arr[I], 2) + ...));
      }(std::make_index_sequence<D>());
    }

    // variadic max base case.
    template<typename T>
    auto vmax(const T arg) {
      return arg;
    }

    // variadic max.
    template<typename T1, typename T2, typename... T>
    auto vmax(const T1 arg1, const T2 arg2, const T... args) {
        return (arg1 > arg2) ? vmax(arg1, args...) : vmax(arg2, args...);
    }

    // L_inf norm
      template<std::size_t D, typename T>
    auto chessboard_distance(const std::array<T, D> arr) {
      return [&]<std::size_t... I>(std::index_sequence<I...>) {
        return vmax(std::abs(arr[I])...);
      }(std::make_index_sequence<D>());
    }

} }






//--------------------------------------------------
// Junkyard


// Eulerian distances (i.e., L_2 norm)
//--------------------------------------------------
//
//template<std::size_t D, typename T = std::enable_if_t<(D==1),int> >
//inline double eulerian_distance(corgi::internals::tuple_of<1, int> rel)
//{
//  return (double) std::abs(std::get<0>(rel));
//}
//
//template<std::size_t D, typename T = std::enable_if_t<(D==2),int> >
//inline double eulerian_distance(corgi::internals::tuple_of<2, int> rel)
//{
//  return sqrt( 
//      (double)std::pow(std::get<0>(rel),2) +
//      (double)std::pow(std::get<1>(rel),2)
//      );
//}
//
//template<std::size_t D, typename T = std::enable_if_t<(D==3),int> >
//inline double eulerian_distance(corgi::internals::tuple_of<3, int> rel)
//{
//  return sqrt( 
//      (double)std::pow(std::get<0>(rel),2) +
//      (double)std::pow(std::get<1>(rel),2) +
//      (double)std::pow(std::get<2>(rel),2)
//      );
//}

