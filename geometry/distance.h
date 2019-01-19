#pragma once

#include <tuple>
#include <cmath>

#include "../internals.h"
#include "utilities.h"


namespace corgi { 
  
  namespace geom {



    // L_1 norm
    //
    // FIXME: C++17 expression not supported by INTEL
    //
    //template<class Tuple>
    //decltype(auto) manhattan_distance(const Tuple& tup)
    //{
    //  auto norm = [](auto const&... e)->decltype(auto) 
    //  { 
    //    return ( std::abs(e) +...); 
    //  };

    //  return static_cast<double>( notstd::apply( norm, tup ) );
    //}

    template<std::size_t D, typename T = std::enable_if_t<(D==1),int> >
    decltype(auto) manhattan_distance(
        const corgi::internals::tuple_of<1, T>& tup)
    {
      return std::abs(std::get<0>(tup));
    }

    template<std::size_t D, typename T = std::enable_if_t<(D==2),int> >
    decltype(auto) manhattan_distance(
        const corgi::internals::tuple_of<2, T>& tup)
    {
      return std::abs(std::get<0>(tup)) + std::abs(std::get<1>(tup));
    }



    // L_2 norm
    //
    // FIXME: C++17 expression not supported by INTEL
    //
    //template<class Tuple>
    //decltype(auto) eulerian_distance(const Tuple& tup)
    //{
    //  auto norm = [](auto const&... e)->decltype(auto) 
    //  { 
    //    return ( (e*e) +...); 
    //  };

    //  return sqrt( static_cast<double>( notstd::apply( norm, tup ) ) );
    //}

    template<std::size_t D, typename T = std::enable_if_t<(D==1),int> >
    decltype(auto) euler_distance(
        const corgi::internals::tuple_of<1, T>& tup)
    {
      return std::abs(std::get<0>(tup));
    }

    template<std::size_t D, typename T = std::enable_if_t<(D==2),int> >
    decltype(auto) euler_distance(
        const corgi::internals::tuple_of<2, T>& tup)
    {
      return std::sqrt(std::pow(std::get<0>(tup),2) + std::pow(std::get<1>(tup),2));
    }



    //template<typename... T>
    //constexpr T abs_max(

    // L_inf norm
    //
    // FIXME: hand coded
    //template<class Tuple>
    //decltype(auto) chessboard_distance(const Tuple& tup)
    //{
    //  auto norm = [](auto const& e, auto const&... args)->decltype(auto) 
    //  { 
    //    return std::max(std::abs(e), args...); 
    //  };

    //  return static_cast<double>( norm(tup) );

    //  //return static_cast<double>( std::max( std::abs( std::initializer_list<tup>() ) ) );
    //}

    template<std::size_t D, typename T = std::enable_if_t<(D==1),int> >
    decltype(auto) chessboard_distance(
        const corgi::internals::tuple_of<1, T>& tup)
    {
      return std::abs(std::get<0>(tup));
    }

    template<std::size_t D, typename T = std::enable_if_t<(D==2),int> >
    decltype(auto) chessboard_distance(
        const corgi::internals::tuple_of<2, T>& tup)
    {
      return std::max( 
          std::abs(std::get<0>(tup)), 
          std::abs(std::get<1>(tup))
          );
    }

    template<std::size_t D, typename T = std::enable_if_t<(D==3),int> >
    decltype(auto) chessboard_distance(
        const corgi::internals::tuple_of<3, T>& tup)
    {
      return 
        std::max(
          std::max( std::abs(std::get<0>(tup)), std::abs(std::get<1>(tup)) ),
          std::abs(std::get<2>(tup))
          );
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

