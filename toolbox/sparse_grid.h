#pragma once

#include <map>
#include <array>
#include <vector>
#include <stdexcept>
#include <type_traits>

#include "../internals.h"


namespace corgi {
  namespace tools {


/// \brief Sparse adaptive grid
//
// indexing is via tuples so the grid can grow to any direction 
// (including -directions) as hash( tuple<xxx> ) is always unique.
//
// Internally data is stored in a map<> resembling a sparse matrix.
//
template< typename T, int D>
class sparse_grid {

  private:

  /// internal data storage
  std::map< corgi::internals::tuple_of<D, size_t>, int> _data;
    
  /// number of elements in each dimension
  std::array<size_t, D> _lengths;


  public:

  /// () referencing
  template<
    typename... Dlen,
    typename = corgi::internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi::internals::are_integral<Dlen...>::value , void> 
  >
  T& operator()(Dlen... indices)
  {
    return _data[ std::tuple<Dlen...>(indices...) ];
  }

  /// const () referencing
  template<
    typename... Dlen,
    typename = corgi::internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi::internals::are_integral<Dlen...>::value , void> 
  >
  const T& operator()(Dlen... indices) const
  {
    return _data[ std::tuple<Dlen...>(indices...) ];
  }


  
  // ctor with grid size
  template<
    typename... Dlen,
    typename = corgi::internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi::internals::are_integral<Dlen...>::value , void> 
  >
  sparse_grid(Dlen... lens) :
    _lengths {{ static_cast<size_t>(lens)... }}
  { }


  // ctor without grid size
  sparse_grid() {};


  /// resize the assumed size
  template< typename... Dlen >
  corgi::internals::enable_if_t<(sizeof...(Dlen) == D) &&
  corgi::internals::are_integral<Dlen...>::value , 
    void> 
  resize(Dlen... _lens) 
  {
    std::array<size_t, D> lens = {{static_cast<size_t>(_lens)...}};
    _lengths = lens;
  }
  

  /// Return object that is contiguous in memory
  std::vector<T> serialize() 
  {
    std::vector<T> ret;

    // TODO: implementation
    // XXX: how to recursively loop over _length?
    //
    //ret.resize(Nx * Ny);
    //for(size_t k=0; k<Nx*Ny; k++) {ret[k] = 0.0;};

    //// get elements from map
    //size_t indx;
    //for(auto const elem: mat) {
    //  indx = Nx*elem.first.second + elem.first.first;
    //  ret[indx] = elem.second;
    //}

    return ret;
  }
  

  template< typename... Dlen >
  corgi::internals::enable_if_t<(sizeof...(Dlen) == D) &&
  corgi::internals::are_integral<Dlen...>::value , 
    void> 
  unpack(std::vector<T>& /*vec*/, Dlen... _lens) 
  {
    std::array<size_t, D> lens = {{static_cast<size_t>(_lens)...}};

    // clear before actually unpacking
    clear();

    // TODO: implementation
    // XXX: how to recursively loop over _length?
    //  size_t k=0;
    //  for(size_t j=0; j<Ny; j++) {
    //    for(size_t i=0; i<Nx; i++) {
    //      if(vec[k] != 0.0) {
    //        mat[ std::make_pair(i,j) ] = vec[k];
    //      }
    //      k++;
    //    }
    //  }

  }


  /// erase the map container
  void clear() {
    _data.clear();
  }


};




  } // end of tools
} // end of corgi
