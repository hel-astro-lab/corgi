#pragma once

#include <vector>
#include <tuple>
#include <cmath>

#include "corgi/internals.h"
#include "corgi/geometry/distance.h"

namespace corgi { namespace ca {





/// Moore neighborhood of different dimensions
/// using type constraints to pick dimensionality specialization

template<std::size_t D>
  requires (D == 1)
std::vector<std::array<int, 1>> moore_neighborhood()
{
  return std::vector<std::array<int, 1>> { { -1 }, { 1 } };
}


template<std::size_t D>
  requires (D == 2)
std::vector<std::array<int, 2> > moore_neighborhood()
{
  std::vector<std::array<int, 2> > ret;
  ret.reserve(8);
  for (int jr=-1; jr<=1; jr++) {
    for (int ir=-1; ir<=1; ir++) {
      if (!( ir == 0 && jr == 0  )) {
        ret.push_back( {ir, jr} );
      }
    }
  }
  return ret;
}

template<std::size_t D>
  requires (D == 3)
std::vector<std::array<int, 3> > moore_neighborhood()
{
  std::vector<std::array<int, 3> > ret;
  ret.reserve(26);
  for (int kr=-1; kr<=1; kr++) {
    for (int jr=-1; jr<=1; jr++) {
      for (int ir=-1; ir<=1; ir++) {
        if (!( ir == 0 && jr == 0 && kr == 0 )) {
          ret.push_back( {ir, jr, kr} );
        }
      }
    }
  }
  return ret;
}


// Box distances
//--------------------------------------------------
template<std::size_t D>
  requires (D == 1)
std::vector<std::array<int, 1> > chessboard_neighborhood(int radius)
{
  std::vector<std::array<int, 1> > ret;
  for (int ir=-radius; ir<=radius; ir++) {
    if (!( ir == 0 )) {
      ret.push_back( {ir} );
    }
  }
  return ret;
}


template<std::size_t D>
  requires (D == 2)
std::vector<std::array<int, 2> > chessboard_neighborhood(int radius)
{
  std::vector<std::array<int, 2> > ret;
  for (int jr=-radius; jr<=radius; jr++) {
    for (int ir=-radius; ir<=radius; ir++)  {
      if (!( ir == 0 && jr == 0 )) {
        ret.push_back( {ir, jr} );
      }
    }
  }
  return ret;
}

template<std::size_t D>
  requires (D == 3)
std::vector<std::array<int, 3> > chessboard_neighborhood(int radius)
{
  std::vector<std::array<int, 3> > ret;
  for (int kr=-radius; kr<=radius; kr++) {
    for (int jr=-radius; jr<=radius; jr++) {
      for (int ir=-radius; ir<=radius; ir++)  {
        if (!( ir == 0 && jr == 0 && kr == 0)) {
          ret.push_back( {ir, jr, kr} );
        }
      }
    }
  }
  return ret;
}

// Spherical distances
//--------------------------------------------------
template<std::size_t D>
  requires (D == 1)
std::vector<std::array<int, 1> > euler_neighborhood(int radius)
{
  return chessboard_neighborhood<1>(radius);
}


template<std::size_t D>
  requires (D == 2)
std::vector<std::array<int, 2> > euler_neighborhood(int radius)
{
  std::vector<std::array<int, 2> > ret;
  for (int jr=-radius; jr<=radius; jr++) {
    for (int ir=-radius; ir<=radius; ir++) {

      auto rel = std::array{ ir,jr };
      if( ::corgi::geom::eulerian_distance<D>(rel) >= (double)radius ) continue;

      if (!( ir == 0 && jr == 0 )) {
        ret.push_back( rel );
      }
    }
  }
  return ret;
}

template<std::size_t D>
  requires (D == 3)
std::vector<std::array<int, 3> > euler_neighborhood(int radius)
{
  std::vector<std::array<int, 3> > ret;
  for (int kr=-radius; kr<=radius; kr++) {
    for (int jr=-radius; jr<=radius; jr++) {
      for (int ir=-radius; ir<=radius; ir++) {

        auto rel = std::array{ ir,jr,kr };
        if( ::corgi::geom::eulerian_distance<D>(rel) >= (double)radius ) continue;

        if (!( ir == 0 && jr == 0 && kr == 0)) {
          ret.push_back( rel );
        }
      }
    }
  }
  return ret;
}




} } // ns corgi::ca
