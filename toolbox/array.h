#pragma once


#include <vector>

namespace toolbox {


template<typename T, size_t N=1>
  class array {

    public:
      typedef T  value_type;
      typedef T* pointer;
      typedef T& reference;

      typedef std::shared_ptr<value_type> shared_pointer;
      typedef array<T,N> self;


    private:

      /// Internal 1D data storage
      std::vector<T> data;

      /// array dimensions
      // std::array<size_t, N> dims;
      std::vector<size_t> dims;


      /// Allocate internal storage array
      void _allocate() {

        size_t ndim = 1;
        for(auto idim : dims) ndim *= idim;

        data.resize(ndim);
      }

    public:

      //--------------------------------------------------
      // constructing/initializing

      /// initialize with variadic args; initialization must occur in ctor
      template <typename ...Args> 
      array(Args ...args) : dims({args...}) {
        static_assert(sizeof...(args) == N, "Incorrect number of arguments");

        _allocate();
      }

      /// initializer list ctor
      array(std::initializer_list<T> args) : dims(args) {
        static_assert(sizeof(args) == N, "Incorrect number of arguments");

        _allocate();
      }


      //--------------------------------------------------
      // accessing/indexing
        
        
      /// clear content
      void clear() {
        data.clear();
      }

      void resize(size_t Ni, ...) {
        va_list vl;
        va_start(vl,Ni);

        // TODO implement allocation
      }


      /// 1D access
      inline reference operator[] (index_t n) const {
        return data[n];
      }

      /// multidimensional access
      reference operator() (index_ptr subs) const {
        return data[ sub2ind(subs, m_size, m_strides) ];
      }

      /// using initializer list, i.e., {1,2,3...}
      reference operator() ( std::initializer_list<index_t> subs ) const {
#ifdef DEBUG
        if (sibs.size() != N) 
          throw std::length_error("Invalid coordinate length");
#endif

        return data[ sub2ind(subs.begin(), m_size, m_strides) ];
      }

      /// coordinate access
      reference operator() (index_t i, ...) const {
        va_list vl;
        va_start(vl,i);

        return data[ (i*m_strides[0]) + sub2ind(vl, m_size, m_strides) ];
      }




  };


}







