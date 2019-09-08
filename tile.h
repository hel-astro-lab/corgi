#pragma once


#include <array>
#include <vector>
#include <tuple>
#include <iostream>
#include <algorithm>

#include "common.h"
#include "internals.h"
#include "cellular_automata.h"

#include <mpi4cpp/mpi.h>


namespace corgi {

namespace mpi = mpi4cpp::mpi;


// Data storage struct for communication members
//
// This is a general D-independent class to make MPI communications
// easier so that we always communicate only this object and use it
// to reconstruct the back.
struct Communication {

  /// my index
  int cid;

  /// (i,j,k) indices
  std::array<int, 3> indices;

  /// MPI rank of who owns me
  int owner;

  /// If I am a virtual tile, who do I share the values the most.
  int top_virtual_owner;

  /// how many times do I have to be sent to others
  int communications = 0;

  /// How many virtual neighbors do I have
  int number_of_virtual_neighbors = 0;

  /// tile boundaries
  std::array<double, 3> mins;
  std::array<double, 3> maxs;
    
  /// tile type listing
  //bool local = false;

  /// my virtual owners; not communicated along the metainfo
  //std::vector<int> virtual_owners;


};



/*! \brief Tile object
 *
 * This is the smallest storage unit of the framework. Internally this should host 
 * a mesh/grid/particles/etc.
 */
template<std::size_t D>
class Tile
{

  public:

    /// unique tile ID
    uint64_t cid;

    // Order here is fixed for mpi_tile_t
    Communication communication;

    // my virtual owners
    std::vector<int> virtual_owners;
    
    /// coarse mpi_grid grid indices
    corgi::internals::tuple_of<D, size_t> index;

    /// Global grid dimensions (needed for wrapping boundaries)
    std::array<size_t, D> lengths;

    /// tile boundaries
    std::array<double, D> mins;
    std::array<double, D> maxs;

      
    // using default ctor
    // TODO: are there side effects?
    Tile() = default;

    /*! \brief *virtual* base class destructor 
     * NOTE: this needs to be virtual so that child classes can be 
     * destroyed.
     */
    virtual ~Tile() = default;

    /// load tile metainfo from Communication object
    void load_metainfo(Communication cm)
    {
      communication = cm;
      cid = cm.cid; 

      // create temporary array of correct size
      std::array<size_t, D> ind2;
      for(size_t i=0; i<D; i++) ind2[i] = static_cast<size_t>(cm.indices[i]);

      // cast into tuple
      index = corgi::internals::into_tuple(ind2);

      for(size_t i=0; i<D; i++) mins[i] = cm.mins[i];
      for(size_t i=0; i<D; i++) maxs[i] = cm.maxs[i];

    }


    /// general N-dim implementation of wrap
    size_t wrap(int ind, size_t d)
    {
      auto N = static_cast<int>(lengths[d]);
      while (ind < 0) { ind += N; }
      while (ind >= N) { ind -= N; }
      return static_cast<size_t>(ind);
    }


    /// return index of tiles in relative to my position
  public:
    template<typename... Indices>
    corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
    corgi::internals::are_integral<Indices...>::value, 
    const corgi::internals::tuple_of<D, size_t> > 
    neighs(Indices... indices_rel)
    {
        std::array<int, D>    rel = {{static_cast<int>(indices_rel)...}};
        //std::array<size_t, D> cur = {{index}};
        auto cur = corgi::internals::into_array(index);

        for(size_t i=0; i<D; i++) {
          cur[i] = static_cast<size_t>(
            wrap( static_cast<int>(rel[i]) + 
                  static_cast<int>(cur[i]), i)
          );
        }

        auto ret = corgi::internals::into_tuple(cur);
        return ret;
    }

    /// auxiliary function to unpack tuples
  private:
    template <size_t... Is>
    const corgi::internals::tuple_of<D, size_t> neighs_impl(
        corgi::internals::tuple_of<D, int>& tuple, 
        std::index_sequence<Is...>)
    {
      return neighs( std::get<Is>(tuple)... );
    }

  public:
    /// unpack tuple into variadic argument list
    template<typename Indices = std::make_index_sequence<D>>
    const corgi::internals::tuple_of<D, size_t> neighs( 
        corgi::internals::tuple_of<D, int>& indices)
    {
        return neighs_impl(indices, Indices{} );
    }

    // end of neighs + auxiliary functions
    //--------------------------------------------------


    /// Return full Moore neighborhood around me
    std::vector< corgi::internals::tuple_of<D, size_t> > nhood()
    {
      std::vector< corgi::internals::tuple_of<D, size_t> > nh;
      for(auto& reli : corgi::ca::moore_neighborhood<D>() ){
        nh.push_back( neighs(reli) );
      }
      return nh;
    }


    // --------------------------------------------------
    // (optional) tile geometry 

    /// set tile minimum limits
    void set_tile_mins(std::array<double, D> bbox)
    {
      mins = std::move(bbox);
      for(size_t i=0; i<D; i++) communication.mins[i] = bbox[i];
    }

    /// set tile minimum limits
    void set_tile_maxs(std::array<double, D> bbox)
    {
      maxs = std::move(bbox);
      for(size_t i=0; i<D; i++) communication.maxs[i] = bbox[i];
    }

    // --------------------------------------------------

    /// send basic information of this (corgi::Tile) to dest
    //mpi::request send(mpi::communicator& comm, int dest)
    //{
    //  return comm.isend(this->communicator, dest);
    //}

    ///// receive basic information of this (corgi::Tile) from orig
    //mpi::request recv(mpi::communicator& comm, int orig)
    //{
    //  return comm.irecv(this->communicator, orig);
    //}

    /// dummy MPI data send function
    virtual std::vector<mpi::request> 
    send_data(
        mpi::communicator& /*comm*/,
        int dest, 
        int /*mode*/,
        int /*tag*/)
    {
      std::vector<mpi::request> reqs;

      std::cout << "send to " << dest << "\n";

      return reqs;
    }


    /// dummy MPI data recv function
    virtual std::vector<mpi::request> 
    recv_data(
        mpi::communicator& /*comm*/,
        int orig, 
        int /*mode*/,
        int /*tag*/)
    {
      std::vector<mpi::request> reqs;

      std::cout << "recv from " << orig << "\n";

      return reqs;
    }

    /// dummy MPI data recv function for extra data 
    virtual std::vector<mpi::request> 
    recv_extra_data(
        mpi::communicator& /*comm*/,
        int orig, 
        int /*mode*/,
        int /*tag*/)
    {
      std::vector<mpi::request> reqs;

      std::cout << "recv from " << orig << "\n";

      return reqs;
    }


    /// Local computational work estimate for this tile
    virtual double get_work()
    {
      return 1.0;
    }



}; // end of Tile class

} // end of namespace corgi
