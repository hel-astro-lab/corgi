#pragma once

#include <map>
#include <vector>
#include <stdexcept>
#include <tuple>


/// Simple sparse matrix class
template <class T>
class SparseGrid{
  public:
    std::map<std::tuple<size_t, size_t>, T> mat;
    size_t Nx;
    size_t Ny;

    /// standard (i,j) syntax
    T& operator()(size_t i, size_t j) { 
      return mat[std::make_tuple(i,j)];
    }

    const T& operator()(size_t i, size_t j) const { 
      return mat[std::make_tuple(i,j)];
    }

    /// std::tuple syntax
    T& operator()(std::tuple<size_t, size_t> indx) { 
      return mat[indx];
    }

    const T& operator()(std::tuple<size_t, size_t> indx) const { 
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

    /// Return object that is contiguous in memory
    std::vector<T> serialize() {

      // create zero array
      std::vector<T> ret;
      ret.resize(Nx * Ny);
      for(size_t k=0; k<Nx*Ny; k++) {ret[k] = 0.0;};

      // get elements from map
      size_t indx;
      for(const auto&& elem: mat) {
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
            mat[ std::make_tuple(i,j) ] = vec[k];
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
