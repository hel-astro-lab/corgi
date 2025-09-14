#pragma once


#include <array>
#include <vector>
#include <tuple>
#include <iostream>
#include <algorithm>
#include <type_traits>

#include "corgi/common.h"
#include "corgi/internals.h"
#include "corgi/cellular_automata.h"

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
//  public std::enable_shared_from_this<Tile<D>>
{

  public:

    /// unique tile ID
    uint64_t cid;

    // Order here is fixed for mpi_tile_t
    Communication communication;

    // my virtual owners
    std::vector<int> virtual_owners;
    
    /// coarse mpi_grid grid indices
    std::array<size_t, D> index;

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
      std::array<std::size_t, D> ind2;
      for(size_t i=0; i<D; i++) ind2[i] = static_cast<std::size_t>(cm.indices[i]);

      index = ind2;

      for(size_t i=0; i<D; i++) mins[i] = cm.mins[i];
      for(size_t i=0; i<D; i++) maxs[i] = cm.maxs[i];

    }

    /// general N-dim implementation of wrap
    size_t wrap(int ind, size_t d) const
    {
      auto N = static_cast<int>(lengths[d]);
      assert(N < 1e6); // FIXME: debug catch preventing ridiculously large while loops
      assert(N > 0);   // FIXME: debug catch preventing ridiculously large while loops

      while (ind < 0) { ind += N; }
      while (ind >= N) { ind -= N; }
      return static_cast<size_t>(ind);
    }

    /// return index of tiles in relative to my position
    template<typename... Indices>
      requires indices_for<D, Indices...>
    std::array<size_t, D> neighs(Indices... indices_rel) const
    {
        std::array<int, D> rel = {{static_cast<int>(indices_rel)...}};
        auto cur = index;

        for(size_t i=0; i<D; i++) {
          cur[i] = static_cast<size_t>(
            wrap(static_cast<int>(rel[i]) + static_cast<int>(cur[i]), i));
        }
        return cur;
    }

    std::array<size_t, D> neighs(const std::array<int, D>& indices) const
    {
        return [&, this]<std::size_t... I>(std::index_sequence<I...>) {
            return neighs(indices[I]...);
        }(std::make_index_sequence<D>());
    }

    // end of neighs + auxiliary functions
    //--------------------------------------------------


    /// Return full Moore neighborhood around me
    std::vector<std::array<std::size_t, D>> nhood()
    {
      std::vector<std::array<std::size_t, D>> nh;
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

    /// dummy pairwise Moore neighborhood communication.
    ///
    /// For each local tile A, corgi::Grid::pairwise_moore_communication
    /// will call A.pairwise_moore_communication(B, dir_to_B, mode)
    /// for each tile B in A's Moore neighborhood. Before this
    /// A.pairwise_moore_communication_prelude(mode) will be called
    /// on each tile and A.pairwise_moore_communication_postlude(mode)
    /// will be called after.
    ///
    /// In addition to local tiles pre- and postlude is called on virtual tiles.
    ///
    /// Order of the calls in corgi::Grid::pairwise_moore_communication
    /// is that Tile::pairwise_moore_communication(B, dir_to_B, mode)
    /// is called for all tiles with same direction before proceeding
    /// to the next direction. There is no guarantees on the order of directions
    /// or the order of tiles.
    virtual void
    pairwise_moore_communication(const Tile& /* other */,
                                 const std::array<int, D> dir_to_other,
                                 const int /* mode */
    ) {
        std::cout << "pairwise moore communication towards ";
        for (const auto x : dir_to_other) { std::cout << x << " "; }
        std::cout << "\n";
    };

    virtual void pairwise_moore_communication_prelude(const int /* mode */) { }
    virtual void pairwise_moore_communication_postlude(const int /* mode */) { }


    /// Local computational work estimate for this tile
    virtual double get_work()
    {
      return 1.0;
    }


}; // end of Tile class

} // end of namespace corgi
