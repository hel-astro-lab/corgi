#pragma once

#include <map>
#include <array>
#include <vector>
#include <stdexcept>
#include <type_traits>

#include "corgi/internals.h"



namespace corgi {
  namespace tools {

namespace internal {

  //template<typename T, int D>
  //void serialize( 
  //  std::map< corgi::internals::tuple_of<D, size_t>, T>& data,
  //  std::vector<T>& ret,
  //  std::array<size_t, 1>& lengths
  //) = delete;
  // 1D

  template<typename T>
  void serialize(
    std::map<std::array<std::size_t, 1>, T>& data,
    std::vector<T>& ret,
    std::array<size_t, 1>& /*lengths*/)
  {
    size_t indx;
    for(auto const& elem: data) {
      indx = std::get<0>(elem.first);
      ret[indx] = elem.second;
    }
  }

  // 2D
  template<typename T>
  void serialize( 
    std::map<std::array<std::size_t, 2>, T>& data,
    std::vector<T>& ret,
    std::array<size_t, 2>& lengths) 
  {
    size_t indx;
    for(auto const& elem: data) {
      indx  =            std::get<0>(elem.first);
      indx += lengths[0]*std::get<1>(elem.first);
      ret[indx] = elem.second;
    }
  }

  // 3D
  template<typename T>
  void serialize( 
    std::map<std::array<std::size_t, 3>, T>& data,
    std::vector<T>& ret,
    std::array<size_t, 3>& lengths) 
  {
    size_t indx;
    for(auto const& elem: data) {
      indx  =                       std::get<0>(elem.first);
      indx +=            lengths[0]*std::get<1>(elem.first);
      indx += lengths[0]*lengths[1]*std::get<2>(elem.first);
      ret[indx] = elem.second;
    }
  }


  // 1D
  template<typename T>
  void deserialize( 
    std::map<std::array<std::size_t, 1>, T>& data,
    std::vector<T>& inc,
    std::array<size_t, 1>& lengths)
  {
    for(size_t i=0; i<lengths[0]; i++) {
      data[ std::array{ i } ] = inc[i];
    }
  }

  // 2D
  template<typename T>
  void deserialize( 
    std::map<std::array<std::size_t, 2>, T>& data,
    std::vector<T>& inc,
    std::array<size_t, 2>& lengths)
  {

    size_t indx;
    for(size_t j=0; j<lengths[1]; j++) {
    for(size_t i=0; i<lengths[0]; i++) {
      indx  = i;
      indx += lengths[0]*j;
      data[ std::array{ i,j } ] = inc[indx];
    }
    }
  }

  // 3D
  template<typename T>
  void deserialize( 
    std::map<std::array<std::size_t, 3>, T>& data,
    std::vector<T>& inc,
    std::array<size_t, 3>& lengths)
  {

    size_t indx;
    for(size_t k=0; k<lengths[2]; k++) {
    for(size_t j=0; j<lengths[1]; j++) {
    for(size_t i=0; i<lengths[0]; i++) {
      indx  = i;
      indx += lengths[0]*j;
      indx += lengths[0]*lengths[1]*k;
      data[ std::array{ i,j,k } ] = inc[indx];
    }
    }
    }
  }
}

//--------------------------------------------------


/// \brief Sparse adaptive grid
//
// indexing is via tuples so the grid can grow to any direction 
// (including -directions) as hash( tuple<xxx> ) is always unique.
//
// Internally data is stored in a map<> resembling a sparse matrix.
//
template< typename T, int D>
class sparse_grid {

  using map_t = std::map<std::array<std::size_t, D>, T>;


  private:

  /// internal data storage
  map_t _data;
    
  /// number of elements in each dimension
  std::array<size_t, D> _lengths;


  public:

  /// () referencing
  template<typename... Dlen>
    requires indices_for<D, Dlen...>
  T& operator()(Dlen... indices)
  {
    return _data[std::array{ indices... }];
  }

  T& operator()(const std::array<std::size_t, D> ind)
  {
    return _data[ ind ];
  }


  /// const () referencing
  template<typename... Dlen>
    requires indices_for<D, Dlen...>
  const T& operator()(Dlen... indices) const
  {
    return _data[ std::array{ indices... } ];
  }

  const T& operator()(std::array<std::size_t, D> ind) const
  {
    return _data[ ind ];
  }


  
  // ctor with grid size
  template<typename... Dlen>
    requires indices_for<D, Dlen...>
  sparse_grid(Dlen... lens) :
    _lengths {{ static_cast<size_t>(lens)... }}
  { }


  // ctor without grid size
  sparse_grid() = default;

  // copy ctor
  sparse_grid(const sparse_grid& p)
  {
    _lengths = p._lengths;
    _data = p._data;
  }

  // move assignment ctor
  sparse_grid& operator= (sparse_grid&& other) noexcept
  {
    //std::move(other._data.begin(), other._data.end(), 
    //    _data.begin());
    //std::copy(other._lengths.begin(), other._lengths.end(),
    //    _lengths.begin());
    _lengths = std::move(other._lengths);
    _data = std::move(other._data);

    return *this;
  }

  virtual ~sparse_grid() = default;


  /// resize the assumed size
  template<typename... Dlen>
    requires indices_for<D, Dlen...>
  void resize(Dlen... _lens)
  {
    std::array<size_t, D> lens = {{static_cast<size_t>(_lens)...}};
    _lengths = lens;
  }

  /// Return object that is contiguous in memory
  std::vector<T> serialize() 
  {
    std::vector<T> ret; 

    int N = 1;
    for (size_t i = 0; i<D; i++) N *= _lengths[i];
    ret.resize(N);

    internal::serialize(_data, ret, _lengths);

    return ret;
  }
  

  //template< typename... Dlen >
  //std::enable_if_t<(sizeof...(Dlen) == D) &&
  //corgi::internals::are_integral<Dlen...>::value , 
  //  void> 
  //unpack(std::vector<T>& /*vec*/, Dlen... _lens) 

  void deserialize(std::vector<T>& vec, std::array<size_t, D> lens)
  {

    // copy new size in
    _lengths = lens;

    int N = 1;
    //std::array<size_t, D> lens = {{static_cast<size_t>(_lens)...}};
    for (size_t i = 0; i<D; i++) N *= _lengths[i];

    // clear before actually unpacking
    clear();

    internal::deserialize(_data, vec, _lengths);
  }


  /// erase the map container
  void clear() {
    _data.clear();
  }


  //-------------------------------------------------- 
  // iterators

  typename map_t::iterator begin()
  { 
    return _data.begin();
  }

  typename map_t::const_iterator cbegin()
  { 
    return _data.cbegin();
  }

  typename map_t::iterator end()
  { 
    return _data.end();
  }

  typename map_t::const_iterator cend()
  { 
    return _data.cend();
  }



};




  } // end of tools
} // end of corgi
