#pragma once

#include <map>
#include <vector>
#include <stdexcept>

#include "../internals.h"


namespace tools {


/// Simple sparse matrix class
//template <class T, std::size_t D>
template < typename T >
class SparseGrid{

  private:
    
    template<typename ... Vals>
    using tuple_pack = typename std::tuple<Vals...>;


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
// try 1

template<typename T, typename... Args>
struct storage{
  std::map<std::tuple<Args...>, T> data;
};

//--------------------------------------------------
// try 2
//
// n-length tuple of type T, i.e., tuple_of<3, int> = tuple<int, int, int>
// see: 
//  - https://stackoverflow.com/questions/38885406/produce-stdtuple-of-same-type-in-compile-time-given-its-length-by-a-template-a
template <size_t I,typename T> 
struct tuple_n{
    template< typename...Args> using type = typename tuple_n<I-1, T>::template type<T, Args...>;
};

template <typename T> 
struct tuple_n<0, T> {
    template<typename...Args> using type = std::tuple<Args...>;   
};
template <size_t I,typename T>  using tuple_of = typename tuple_n<I,T>::template type<>;



//--------------------------------------------------
template< typename T, int D>
class sparse_grid {

  /*
  template<
    typename... Dlen,
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  using tuple_map = std::map< std::tuple<Dlen...>, T>;
  */

  // works but does not accomplish anything
  //template<typename... Args>
  //using tuple_map = std::map<std::tuple<Args...>, T>;

  //---

  template<int Dim, typename... Args>
  using check_index_length_t = 
    typename  corgi_internals::enable_if_t<(sizeof...(Args) == Dim) &&
              corgi_internals::are_integral<Args...>::value, void>::type;

  template<typename... Dlen, check_index_length_t<D, Dlen...>>
  storage<T, Dlen...> data2();

  //---

  /*
  template<
    typename... Dlen,
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  using tuple_map2 = std::map<std::tuple<Dlen...>, T>;
  */



  template<
    typename... Dlen,
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  storage<T, Dlen...> data();


  //--------------------------------------------------
  // try 2

  std::map<tuple_of<D, size_t>, int> data3;


  //--------------------------------------------------
  public:

  template<
    typename... Dlen,
    //check_index_length_t<D, Dlen...>
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  //T& operator()(std::tuple<Dlen...>& indices)
  T& operator()(Dlen... indices)
  {
    //return (*this->data2).data[indices];
    //return *(this->data2).data[ std::tuple<Dlen...>(indices...) ];
    //return *(this->data2).data->operator[]( std::tuple<Dlen...>(indices...) );
    return data3[ std::tuple<Dlen...>(indices...) ];
  }

  template<
    typename... Dlen,
    //check_index_length_t<D, Dlen...>
    typename = corgi_internals::enable_if_t<(sizeof...(Dlen) == D) &&
               corgi_internals::are_integral<Dlen...>::value , void> 
  >
  //const T& operator()(std::tuple<Dlen...>& indices) const
  const T& operator()(Dlen... indices) const
  {
    //return (*this->data2).data[indices];
    //return *(this->data2).data->operator[]( std::tuple<Dlen...>(indices...) );
    return data3[ std::tuple<Dlen...>(indices...) ];
  }


};




}
