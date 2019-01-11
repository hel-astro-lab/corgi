#pragma once

#include <vector>
#include <tuple>
#include <cmath>

#include "internals.h"
#include "geometry/distance.h"

namespace corgi { namespace ca {





/// Moore neighborhood of different dimensions
// using SFINAE to pick dimensionality specialization

template<std::size_t D, typename T = std::enable_if_t<(D==1),int> >
std::vector< corgi::internals::tuple_of<D, T> > moore_neighborhood()
{
  std::vector<  corgi::internals::tuple_of<D,int> > ret = {
    {std::make_tuple(-1), std::make_tuple(+1)}
  };
  return ret;
}


template<std::size_t D, typename T = std::enable_if_t<(D==2),int> >
std::vector< corgi::internals::tuple_of<2,int> > moore_neighborhood()
{
  //std::vector<  corgi::internals::tuple_of<2,int> > ret = {
  //  { std::make_tuple(-1, -1), 
  //    std::make_tuple(-1,  0), 
  //    std::make_tuple(-1, +1), 
  //    std::make_tuple( 0, +1), 
  //    std::make_tuple(+1, +1), 
  //    std::make_tuple( 0, +1), 
  //    std::make_tuple(+1, -1), 
  //    std::make_tuple( 0, -1) }
  //};
  std::vector<  corgi::internals::tuple_of<2,int> > ret;
  ret.reserve(8);
  for (int jr=-1; jr<=1; jr++) {
    for (int ir=-1; ir<=1; ir++) {
      if (!( ir == 0 && jr == 0  )) {
        ret.push_back( std::make_tuple(ir, jr) );
      }
    }
  }
  return ret;
}

template<std::size_t D, typename T = std::enable_if_t<(D==3),int> >
std::vector< corgi::internals::tuple_of<3,int> > moore_neighborhood()
{
  std::vector<  corgi::internals::tuple_of<3,int> > ret;
  ret.reserve(26);
  for (int kr=-1; kr<=1; kr++) {
    for (int jr=-1; jr<=1; jr++) {
      for (int ir=-1; ir<=1; ir++) {
        if (!( ir == 0 && jr == 0 && kr == 0 )) {
          ret.push_back( std::make_tuple(ir, jr, kr) );
        }
      }
    }
  }
  return ret;
}


// Box distances
//--------------------------------------------------
template<std::size_t D, typename T = std::enable_if_t<(D==1),int> >
std::vector< corgi::internals::tuple_of<1,int> > chessboard_neighborhood(int radius)
{
  std::vector<  corgi::internals::tuple_of<1,int> > ret;
  for (int ir=-radius; ir<=radius; ir++) {
    if (!( ir == 0 )) {
      ret.push_back( std::make_tuple(ir) );
    }
  }
  return ret;
}


template<std::size_t D, typename T = std::enable_if_t<(D==2),int> >
std::vector< corgi::internals::tuple_of<2,int> > chessboard_neighborhood(int radius)
{
  std::vector<  corgi::internals::tuple_of<2,int> > ret;
  for (int jr=-radius; jr<=radius; jr++) {
    for (int ir=-radius; ir<=radius; ir++)  {
      if (!( ir == 0 && jr == 0 )) {
        ret.push_back( std::make_tuple(ir, jr) );
      }
    }
  }
  return ret;
}


// Spherical distances
//--------------------------------------------------
template<std::size_t D, typename T = std::enable_if_t<(D==1),int> >
std::vector< corgi::internals::tuple_of<1,int> > euler_neighborhood(int radius)
{
  return chessboard_neighborhood<1>(radius);
}


template<std::size_t D, typename T = std::enable_if_t<(D==2),int> >
std::vector< corgi::internals::tuple_of<2,int> > euler_neighborhood(int radius)
{
  std::vector<  corgi::internals::tuple_of<2,int> > ret;
  for (int jr=-radius; jr<=radius; jr++) {
    for (int ir=-radius; ir<=radius; ir++) {

      auto rel = std::make_tuple(ir,jr);
      if( ::corgi::geom::eulerian_distance(rel) >= (double)radius ) continue;

      if (!( ir == 0 && jr == 0 )) {
        ret.push_back( rel );
      }
    }
  }
  return ret;
}





} } // ns corgi::ca
