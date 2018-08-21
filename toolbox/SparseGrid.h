#pragma once

#include <map>
#include <array>
#include <vector>
#include <stdexcept>
#include <type_traits>

#include "../internals.h"


namespace tools {


/// Simple sparse matrix class
//template <class T, std::size_t D>
template < typename T >
class SparseGrid{


  public:
    std::map<std::pair<size_t, size_t>, T> mat;

    size_t Nx;
    size_t Ny;

    /// standard (i,j) syntax
    T& operator()(size_t i, size_t j) { 
      return mat[std::make_pair(i,j)];
    }
    const T& operator()(size_t i, size_t j) const { 
      return mat[std::make_pair(i,j)];
    }

    /// std::pair syntax
    T& operator()(std::pair<size_t, size_t> indx) { 
      return mat[indx];
    }
    const T& operator()(std::pair<size_t, size_t> indx) const { 
      return mat[indx];
    }

    /// Default constructor with no pre-assumed size
    SparseGrid() {
      Nx = 0;
      Ny = 0;
    };

    /*! initialize with information about *assumed* size
     * NOTE: size can be adapted dynamically later on though
     */
    SparseGrid(size_t nx, size_t ny) {
      Nx = nx;
      Ny = ny;
    }


    /// resize existing grid
    void resize(size_t nx, size_t ny) {
      Nx = nx;
      Ny = ny;
    }

    //XXX ugly hack to make this work, for now
    //void resize(std::initializer_list<size_t> lens) {
    //  auto it=lens.begin();
    //  Nx = *it;

    //  it++;
    //  Ny = *it;
    //}
    

    /// Return object that is contiguous in memory
    std::vector<T> serialize() {

      // create zero array
      std::vector<T> ret;
      ret.resize(Nx * Ny);
      for(size_t k=0; k<Nx*Ny; k++) {ret[k] = 0.0;};

      // get elements from map
      size_t indx;
      for(auto const elem: mat) {
        indx = Nx*elem.first.second + elem.first.first;
        ret[indx] = elem.second;
      }

      return ret;
    }

    /// Unpack dense vector into sparse grid
    void unpack(std::vector<T> vec, size_t Nx, size_t Ny) {
      if (vec.size() != Nx*Ny) {
        throw std::range_error("size of vector does not match given input");
      }

      // clear mat before actually unpacking
      clear();

      size_t k=0;
      for(size_t j=0; j<Ny; j++) {
        for(size_t i=0; i<Nx; i++) {
          if(vec[k] != 0.0) {
            mat[ std::make_pair(i,j) ] = vec[k];
          }
          k++;
        }
      }

    }

    /// erase the map container
    void clear() {
      mat.clear();
    }

};




//--------------------------------------------------
template< typename T, int D>
class sparse_grid {

  private:

  /// internal data storage
  std::map< corgi_internals::tuple_of<D, size_t>, int> _data;
    
  /// number of elements in each dimension
  std::array<size_t, D> _lengths;



  public:

  /// () referencing
  template<
    typename... Dlen,
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  T& operator()(Dlen... indices)
  {
    return _data[ std::tuple<Dlen...>(indices...) ];
  }

  /// const () referencing
  template<
    typename... Dlen,
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  const T& operator()(Dlen... indices) const
  {
    return _data[ std::tuple<Dlen...>(indices...) ];
  }


  
  // ctor with grid size
  template<
    typename... Dlen,
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  sparse_grid(Dlen... lens) :
    _lengths {{ static_cast<size_t>(lens)... }}
  { }


  // ctor without grid size
  sparse_grid() {};


  /// resize the assumed size
  template< typename... Dlen >
  corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
  corgi_internals::are_integral<Dlen...>::value , 
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

    return ret;
  }
  

  template< typename... Dlen >
  corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
  corgi_internals::are_integral<Dlen...>::value , 
    void> 
  unpack(std::vector<T>& /*vec*/, Dlen... _lens) 
  {
    std::array<size_t, D> lens = {{static_cast<size_t>(_lens)...}};

    // clear before actually unpacking
    clear();

    // TODO: implementation
    // XXX: how to recursively loop over _length?

  }


  /// erase the map container
  void clear() {
    _data.clear();
  }


};




}
