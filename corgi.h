#pragma once

// #include <fmt/format.h>
// #include <fmt/format.cc>
// #include <fmt/string.h>
// #include <fmt/ostream.h>

#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <cassert>
#include <initializer_list>
#include <sstream>
#include <utility>

#include "internals.h"
#include "toolbox/sparse_grid.h"
#include "tile.h"

//#include "mpi.h"
#include <mpi4cpp/mpi.h>
#include "communication.h"


namespace corgi {

namespace mpi = mpi4cpp::mpi;


/*! Individual node object that stores patches of grid in it.
 *
 * See:
 * - https://github.com/maddouri/hyper_array/
 * - https://github.com/astrobiology/orca_array
*/

template<std::size_t D>
class Node
{


  public:
      
  // --------------------------------------------------
  // definitions
  using size_type  = std::size_t;
  using index_type = std::size_t;
  using float_type = double;


  protected:

  // --------------------------------------------------
  /// number of elements in each dimension
  ::std::array<size_type, D> _lengths;

  /// start coordinates of each dimension
  ::std::array<float_type, D> _mins;
    
  /// ending coordinates of each dimension
  ::std::array<float_type, D> _maxs;

  /*! Global large scale block grid where information
   * of all the mpi processes are stored
   */
  corgi::tools::sparse_grid<int, D> _mpiGrid;

  // --------------------------------------------------
  private:
    
  // Mappings
  using TileID_t = uint64_t;
  using Tile_t   = corgi::Tile<D>;
  using TilePtr  = std::shared_ptr<Tile_t>;
  using TileMap  = std::unordered_map<TileID_t, TilePtr>;



  public:

  /// Map with tileID & tile data
  TileMap tiles;

  // XXX: FIXME: compatibility
  int rank = 0;


  public:
  // --------------------------------------------------
  // Python bindings for mpiGrid

  // get element
  template<typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value, int > 
  pyGetMpiGrid(Indices... indices)  /*const*/
  {
    return _mpiGrid(indices...);
  }


  // set element
  template<typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value, void > 
  pySetMpiGrid(int val, Indices... indices) {
    _mpiGrid(indices...) = val;
  }



  public:

  // --------------------------------------------------
  // constructors

  /// mpi environment
  mpi::environment env;

  /// mpi communicator
  mpi::communicator comm;
    
  /// Uninitialized dimension lengths
  Node() :
    env(),
    comm()
  {};

  /// copy-constructor
  //Node(const Node& /*other*/) {};
  
  /// move constructor
  //Node(Node&& /*other*/) {}; // use std::move()
   
  /// set dimensions during construction time
  template<
    typename... DimensionLength,
    typename = corgi::internals::enable_if_t< (sizeof...(DimensionLength) == D) && 
               corgi::internals::are_integral<DimensionLength...>::value, void
    >
  > 
  Node(DimensionLength... dimensionLengths) :
    _lengths {{static_cast<size_type>(dimensionLengths)...}},
    _mpiGrid(dimensionLengths...),
    env(),
    comm()
  { }


  /*
   * try specializing handy shortcuts to symmetrize construction always assuming 3D input
  template< typename = corgi::internals::enable_if_t< (D == 1), void > > 
  Node(size_t i, size_t j, size_t k) :
    _lengths {{i}},
    _mpiGrid({{i}})
  { }

  template< typename = corgi::internals::enable_if_t< (D == 2), void > > 
  Node(size_t i, size_t j, size_t k) :
    _lengths {{i, j}},
    _mpiGrid({{i, j}})
  { }
  */
  

  /// Deallocate and free everything
  virtual ~Node() = default;

  public:
  
  // --------------------------------------------------
  // assignments
    
  /// copy assignment
  //Node& operator=(const Node& other)
  //{
  //  _lengths = other._lengths;

  //  return *this;
  //}


  /// move assignment
  //Node& operator=(const Node&& other)
  //{
  //  _lengths   = std::move(other._lengths);

  //  return *this;
  //}




  //public:
  // --------------------------------------------------
  // iterators
  /*
          iterator         begin()         noexcept { return iterator(data());                }
    const_iterator         begin()   const noexcept { return const_iterator(data());          }
          iterator         end()           noexcept { return iterator(data() + size());       }
    const_iterator         end()     const noexcept { return const_iterator(data() + size()); }
          reverse_iterator rbegin()        noexcept { return reverse_iterator(end());         }
    const_reverse_iterator rbegin()  const noexcept { return const_reverse_iterator(end());   }
          reverse_iterator rend()          noexcept { return reverse_iterator(begin());       }
    const_reverse_iterator rend()    const noexcept { return const_reverse_iterator(begin()); }
    const_iterator         cbegin()  const noexcept { return const_iterator(data());          }
    const_iterator         cend()    const noexcept { return const_iterator(data() + size()); }
    const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(end());   }
    const_reverse_iterator crend()   const noexcept { return const_reverse_iterator(begin()); }
  */



  public:
  // --------------------------------------------------
  // access grid configuration

  /// number of dimension
  static constexpr size_type dims() noexcept { return D; }

  /// length (run-time)
  size_type len(const size_type i) const
  {
    assert(i < D);
    return _lengths[i];
  }

  /// reference to the _lengths array
  const ::std::array<size_type, D>& lens() const noexcept
  {
    return _lengths;
  }
  
  /// starting location of i:th dimension (run-time)
  float_type min(const size_type i) const
  {
    assert(i < D);
    return _mins[i];
  }
    
  /// ending location of i:th dimension (run-time)
  float_type max(const size_type i) const
  {
    assert(i < D);
    return _maxs[i];
  }

  /// reference to the _mins array
  const ::std::array<float_type, D>& mins() const noexcept
  {
    return _mins;
  }
    
  /// reference to the _maxs array
  const ::std::array<float_type, D>& maxs() const noexcept
  {
    return _maxs;
  }


  // --------------------------------------------------
  // indexing
  private:
  
  template <typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value,
        ::std::array<index_type, D>>
    _validate_index_range(Indices... indices) const
  {
    ::std::array<index_type, D> index_array = {{static_cast<index_type>(indices)...}};

    // check all indices and prepare an exhaustive report (in oss)
    // if some of them are out of bounds
    std::ostringstream oss;
    for (index_type i = 0; i < D; ++i)
    {
      if ((index_array[i] >= _lengths[i]) || (index_array[i] < 0))
      {
        oss << "Index #" << i << " [== " << index_array[i] << "]"
          << " is out of the [0, " << (_lengths[i]-1) << "] range. ";
      }
    }

    // if nothing has been written to oss then all indices are valid
    assert(oss.str().empty());
    return index_array;
  }

  /*! Computes the index coefficients assuming column-major order
   *
   *  what we compute:
   *        \f[
   *            \begin{cases}
   *            C_i = \prod_{j=i+1}^{n-1} L_j
   *            \\
   *            \begin{cases}
   *                i   &\in [0, \text{Dimensions - 1}] \\
   *                C_i &: \text{\_coeffs[i]}           \\
   *                L_j &: \text{\_lengths[j]}
   *            \end{cases}
   *            \end{cases}
   *        \f]
   *
   *  For row-major switch to:
   *  coeffs[i] = ct_accumulate(dimensionLengths, i + 1, Dimensions - i - 1,
   *                                        static_cast<size_type>(1),
   *                                        ct_prod<size_type>);
   *
   */
  std::array<size_type, D>
  compute_index_coeffs(const ::std::array<size_type, D>& dimensionLengths) const noexcept
  {
      std::array<size_type, D> coeffs;
      for (size_type i = 0; i < D; ++i)
      {
          coeffs[i] = corgi::internals::ct_accumulate(
              dimensionLengths,
              0,
              i,
              static_cast<size_type>(1),
              corgi::internals::ct_prod<size_type>);
      }
      return coeffs;
  }


  /// Actual Morton Z-ordering from index list
  // 
  // what we compute: coeff . indices
  //
  // i.e., inner product of accumulated coefficients vector and index vector
  constexpr index_type 
  _compute_index(
      const ::std::array<index_type, D>& index_array) const noexcept
  {
    return corgi::internals::ct_inner_product(
        compute_index_coeffs(_lengths), 0,
        index_array, 0, D,
        static_cast<index_type>(0),
        corgi::internals::ct_plus<index_type>,
        corgi::internals::ct_prod<index_type>);
  }

  public:
    
  /// tile IDs
  template<typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value, index_type > 
  id(Indices... indices) const
  {
    return _compute_index( _validate_index_range(indices...) );
  }


  /// auxiliary function to unpack tuples
  template <size_t... Is>
  index_type id_impl(
      corgi::internals::tuple_of<D, size_t>& tuple, 
      std::index_sequence<Is...>)
  {
    return id( std::get<Is>(tuple)... );
  }

  /// unpack tuple into variadic argument list
  template<typename Indices = std::make_index_sequence<D>>
  index_type id( corgi::internals::tuple_of<D, size_t>& indices)
  {
      return id_impl(indices, Indices{} );
  }
  

  public:

  // --------------------------------------------------
  // apply SFINAE to create some shortcuts (when appropriate) 
  // NOTE: valid up to D=3 with x/y/z

  // return global grid sizes
  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=1), T> 
  getNx() { return _lengths[0]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=2), T> 
  getNy() { return _lengths[1]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=3), T> 
  getNz() { return _lengths[2]; }


  // return global grid limits
  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=1), T> 
  getXmin() { return _mins[0]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=2), T> 
  getYmin() { return _mins[1]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=3), T> 
  getZmin() { return _mins[2]; }


  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=1), T> 
  getXmax() { return _maxs[0]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=2), T> 
  getYmax() { return _maxs[1]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=3), T> 
  getZmax() { return _maxs[2]; }


  /// Set physical grid size
  void setGridLims(
      const ::std::array<float_type, D>& mins,
      const ::std::array<float_type, D>& maxs
      )
  {
    //_mins = std::move(mins);
    //_maxs = std::move(maxs);

    // explicitly avoid move semantics and copy
    // this is to make sure that python garbage collector
    // does not mess things up
    _mins = mins;
    _maxs = maxs;
  }



  public:
  // --------------------------------------------------
  // Tile addition etc. manipulation
    
  /// Add local tile to the node
  // void addTile(Tile& tile) {
  void addTile(
    TilePtr tileptr,
    corgi::internals::tuple_of<D, size_t> indices
    )
  {

    // claim unique ownership of the tile (for unique_ptr)
    // std::unique_ptr<corgi::Tile> tileptr = std::make_unique<corgi:Tile>(tile);
    // TilePtr tileptr = std::make_unique<Tile_t>(tile);
    
    // calculate unique global tile ID
    uint64_t cid = id( indices );

    // Erase any existing tiles to avoid emplace of doing nothing TODO: is this correct?
    tiles.erase(cid);

    tileptr->index               = indices;
    tileptr->cid                 = cid;
    tileptr->communication.cid   = cid;
    tileptr->communication.owner = comm.rank();
    tileptr->communication.local = true; //TODO Catch error if tile is not already mine?
    tileptr->lengths             = _lengths;


    // copy indices from tuple into D=3 array in Communication obj
    auto tmp = corgi::internals::into_array(indices);
    for(size_t i=0; i<D; i++) tileptr->communication.indices[i] = tmp[i];

    // tiles.emplace(cid, std::move(tileptr)); // unique_ptr needs to be moved
    tiles.emplace(cid, tileptr); // NOTE using c++14 emplace to avoid copying
    //tiles.insert( std::make_pair(cid, tileptr) ); // NOTE using c++14 emplace to avoid copying
    //tiles[cid] = tileptr;
    //tiles[cid] = tileptr;
      
    // add to my internal listing
    _mpiGrid( indices ) = comm.rank();
  }


  /// Shortcut for creating raw tiles with only the internal meta info.
  // to be used with message passing (w.r.t addTile that is for use with initialization)
  void createTile(Communication& cm)
  {
    auto tileptr = std::make_shared<Tile_t>();
    tileptr->load_metainfo(cm);

    // additional node info
    tileptr->lengths = _lengths;
    // owner
    // local

    // add
    tiles.emplace(cm.cid, tileptr); // NOTE using c++14 emplace to avoid copying
    _mpiGrid( tileptr->index ) = cm.owner;
  }



  /*! Return a vector of tile indices that fulfill a given criteria.  */
  std::vector<uint64_t> getTileIds(
      const bool sorted=false ) {
    std::vector<uint64_t> ret;

    for (auto& it: tiles) ret.push_back( it.first );

    // optional sort based on the tile id
    if (sorted && !ret.empty()) {
      std::sort(ret.begin(), ret.end());
    }

    return ret;
  }

  /*! \brief Get individual tile (as a reference)
   *
   * NOTE: from StackOverflow (recommended getter method):
   * OtherType& get_othertype(const std::string& name)
   * {
   *     auto it = otMap.find(name);
   *     if (it == otMap.end()) throw std::invalid_argument("entry not found");
   *     return *(it->second);
   * }
   *
   * This way map retains its ownership of the tile and we avoid giving pointers
   * away from the Class.
   */
  Tile_t& getTile(const uint64_t cid) {
    auto it = tiles.find(cid);
    if (it == tiles.end()) { throw std::invalid_argument("tile entry not found"); }

    return *(it->second);
  }

  template<typename... Indices>
    corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
    corgi::internals::are_integral<Indices...>::value, 
  Tile_t&>
  getTileInd(const Indices... indices) {
    uint64_t cid = id(indices...);
    return getTile(cid);
  }

  /// \brief Get individual tile (as a pointer)
  TilePtr getTilePtr(const uint64_t cid) {
    auto it = tiles.find(cid);
    if (it == tiles.end()) { throw std::invalid_argument("entry not found"); }
    return it->second;
  }

  template<typename... Indices>
    corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
    corgi::internals::are_integral<Indices...>::value, 
  TilePtr>
  getTilePtrInd(const Indices... indices)
  {
    uint64_t cid = id(indices...);
    return getTilePtr(cid);
  }

  TilePtr getTilePtr(const std::tuple<size_t> ind) {
    size_t i = std::get<0>(ind);
    return getTilePtrInd(i);
  }

  TilePtr getTilePtr(const std::tuple<size_t, size_t> ind) {
    size_t i = std::get<0>(ind);
    size_t j = std::get<1>(ind);
    return getTilePtrInd(i, j);
  }

  TilePtr getTilePtr(const std::tuple<size_t, size_t, size_t> ind) {
    size_t i = std::get<0>(ind);
    size_t j = std::get<1>(ind);
    size_t k = std::get<2>(ind);
    return getTilePtrInd(i, j, k);
  }


  /// Return all local tiles
  std::vector<uint64_t> getLocalTiles(
      const bool sorted=false ) {

    std::vector<uint64_t> tile_list = getTileIds(sorted);

    size_t i = 0, len = tile_list.size();
    while (i < len) {
      if (!tiles.at( tile_list[i] )->communication.local) {
        std::swap(tile_list[i], tile_list.back());
        tile_list.pop_back();
        len -= 1;
      } else {
        i++;
      }
    }

    return tile_list;
  }


  /// Return all tiles that are of VIRTUAL type.
  std::vector<uint64_t> getVirtuals(
      const bool sorted=false ) {
    std::vector<uint64_t> tile_list = getTileIds(sorted);

    size_t i = 0, len = tile_list.size();
    while (i < len) {
      if (tiles.at( tile_list[i] )->communication.local) {
        std::swap(tile_list[i], tile_list.back());
        tile_list.pop_back();
        len -= 1;
      } else {
        i++;
      }
    }

    return tile_list;
  }


  // /// Check if we have a tile with the given index
  bool isLocal(uint64_t cid) {
    bool local = false;

    // Do we have it on the list=
    if (tiles.count( cid ) > 0) {
      // is it local (i.e., not virtual)
      if ( tiles.at(cid)->communication.local ) {
        local = true;
      }
    }

    return local;
  }

  // // TODO: relative indexing w.r.t. given tile
  // // std::tuple<size_t, size_t> get_neighbor_index(corgi::Tile, int i, int j) {
  // //     return c.neighs( std::make_tuple(i,j) );
  // // }
  // // TODO: get_neighbor_tile(c, i, j)


  // std::vector<int> virtualNeighborhood(uint64_t cid) {

  //   auto c = getTileData(cid);
  //   std::vector< std::tuple<size_t, size_t> > neigs = c->nhood();
  //   std::vector<int> virtual_owners;
  //   for (auto indx: neigs) {

  //     /* TODO: check boundary tiles here; 
  //      * now we assume periodicity in x and y
  //      if (std::get<0>(indx) == ERROR_INDEX ||
  //      std::get<1>(indx) == ERROR_INDEX) {
  //      continue;
  //      }
  //      */

  //     // Get tile id from index notation
  //     size_t i = std::get<0>(indx);
  //     size_t j = std::get<1>(indx);
  //     uint64_t cid = id(i, j);

  //     if (!isLocal( cid )) {
  //       int whoami = _mpiGrid(indx); 
  //       virtual_owners.push_back( whoami );
  //     }
  //   }

  //   return virtual_owners;
  // }


  // // Number of virtual neighbors that the tile might have.
  // /*
  //    size_t number_of_virtual_neighborhood(corgi::Tile c) {
  //    return virtual_neighborhood(c).size();
  //    }
  //    */


  // /*! Analyze my local boundary tiles that will be later on
  //  * send to the neighbors as virtual tiles. 
  //  *
  //  * This is where the magic happens and we analyze what and who to send to.
  //  * These values *must* be same for everybody, this is why we use
  //  * mode of the owner list and in case of conflict pick the smaller value.
  //  * This way everybody knows what to expect and we avoid creating conflicts 
  //  * in communication. This information is then being sent to other processes 
  //  * together with the tiles and is analyzed there by others inside the
  //  * `rank_virtuals` function.
  //  * */
  // void analyzeBoundaryTiles() {

  //   for (auto cid: getTiles()) {
  //     std::vector<int> virtual_owners = virtualNeighborhood(cid);
  //     size_t N = virtual_owners.size();

  //     // If N > 0 then this is a boundary tile.
  //     // other criteria could also apply but here we assume
  //     // neighborhood according to spatial distance.
  //     if (N > 0) {

  //       /* Now we analyze `owner` vector as:
  //        * - sort the vector
  //        * - compute mode of the list to see who owns most of the
  //        * - remove repeating elements creating a unique list. */

  //       // sort
  //       std::sort( virtual_owners.begin(), virtual_owners.end() );

  //       // compute mode by creating a frequency array
  //       // NOTE: in case of same frequency we implicitly pick smaller rank
  //       int max=0, top_owner = virtual_owners[0];
  //       for(size_t i=0; i<virtual_owners.size(); i++) {
  //         int co = (int)count(virtual_owners.begin(), 
  //             virtual_owners.end(), 
  //             virtual_owners[i]);
  //         if(co > max) {      
  //           max = co;
  //           top_owner = virtual_owners[i];
  //         }
  //       } 

  //       // remove duplicates
  //       virtual_owners.erase( unique( virtual_owners.begin(), 
  //             virtual_owners.end() 
  //             ), virtual_owners.end() );


  //       // update tile values
  //       auto c = getTileData(cid);
  //       c->top_virtual_owner = top_owner;
  //       c->communications    = virtual_owners.size();
  //       c->number_of_virtual_neighbors = N;

  //       if (std::find( send_queue.begin(),
  //             send_queue.end(),
  //             cid) == send_queue.end()
  //          ) {
  //         send_queue.push_back( cid );
  //         send_queue_address.push_back( virtual_owners );
  //       }
  //     }
  //   }
  // }


  // /// Clear send queue, issue this only after the send has been successfully done
  void clear_send_queue() {
    send_queue.clear();
    send_queue_address.clear();
  }



  // // --------------------------------------------------
  // // Send queues etc.
  //   
  // /// list of tile id's that are to be sent to others
  std::vector<uint64_t> send_queue;

  // /// list containing lists to where the aforementioned send_queue tiles are to be sent
  std::vector< std::vector<int> > send_queue_address;


  // public:
  // // -------------------------------------------------- 

  std::vector<MPI_Request> sent_info_messages;
  std::vector<MPI_Request> sent_tile_messages;

  std::vector<MPI_Request> recv_info_messages;
  std::vector<MPI_Request> recv_tile_messages;


  // /// Broadcast master ranks mpiGrid to everybody
  void bcastMpiGrid() {

    // total size
    int N = 1;
    for (size_t i = 0; i<D; i++) N *= _lengths[i];
    std::vector<int> tmp;

    if (comm.rank() == 0) {
      tmp = _mpiGrid.serialize();
    } else {
      tmp.resize(N);
      for(int k=0; k<N; k++) {tmp[k] = -1.0;};
    }

    MPI_Bcast(&tmp[0],
        N, 
        MPI_INT, 
        0, 
        MPI_COMM_WORLD
        );

    // unpack
    if(comm.rank() != 0) {
      _mpiGrid.deserialize(tmp, _lengths);
    }
  }

  // /// Issue isends to everywhere
  // // First we send a warning message of how many tiles to expect.
  // // Based on this the receiving side can prepare accordingly.
  // void communicateSendTiles() {

  //   sent_info_messages.clear();
  //   sent_tile_messages.clear();
  //   int j = 0;

  //   for (int dest = 0; dest<Nrank; dest++) {
  //     if(dest == rank) { continue; } // do not send to myself

  //     int i = 0;
  //     std::vector<int> to_be_sent;
  //     for (std::vector<int> address: send_queue_address) {
  //       if( std::find( address.begin(),
  //             address.end(),
  //             dest) != address.end()) {
  //         to_be_sent.push_back( i );
  //       }
  //       i++;
  //     }

  //     // initial message informing how many tiles are coming
  //     // TODO: this whole thing could be avoided by using 
  //     // MPI_Iprobe in the receiving end. Maybe.
  //     uint64_t Nincoming_tiles = uint64_t(to_be_sent.size());

  //     MPI_Request req;
  //     sent_info_messages.push_back( req );

  //     MPI_Isend(
  //         &Nincoming_tiles, 
  //         1,
  //         MPI_UNSIGNED_LONG_LONG,
  //         dest,
  //         commType::NCELLS,
  //         comm,
  //         &sent_info_messages[j] 
  //         );
  //     j++;
  //   }


  //   // send the real tile data now
  //   // We optimize this by only packing the tile data
  //   // once, and then sending the same thing to everybody who needs it.
  //   int i = 0;
  //   for (auto cid: send_queue) {
  //     sendTileData( cid, send_queue_address[i] );
  //     i++;
  //   }

  // }

  /// Send individual tile to dest
  // NOTE: we bounce sending back to tile members,
  //       this way they can be extended for different types of send.
  void send_tile(uint64_t cid, int dest)
  {
    mpi::request req;

    auto& tile = getTile(cid);
    //std::cout << comm.rank() << ": sending cid" << cid << "/" << tile.communication.cid << "\n";
    req = comm.isend(dest, 0, tile.communication);

    // FIXME and make non-blocking
    req.wait();

    return;
  }

  void recv_tile(int orig)
  {
    mpi::request req;

    Communication rcom;
    req = comm.irecv(orig, 0, rcom);

    // FIXME and make non-blocking
    req.wait();

    //std::cout << comm.rank() << ":"
    //  << "cid: " << rcom.cid
    //  << "ind: " << rcom.indices[0] << " " << rcom.indices[1] << " " << rcom.indices[2]
    //  << "owner: " << rcom.owner
    //  << "topo: " << rcom.top_virtual_owner
    //  << "comms: " << rcom.communications
    //  << "numv: " << rcom.number_of_virtual_neighbors
    //  << "mins: " << rcom.mins[0] << " " << rcom.mins[1] << " " << rcom.mins[2]
    //  << "maxs: " << rcom.maxs[0] << " " << rcom.maxs[1] << " " << rcom.maxs[2]
    //  << "\n";

    // next need to build tile
    rcom.local = false; // received tiles are automatically virtuals
    createTile(rcom);

    return;
  }


  // /// Pack tile and send to everybody on the dests list
  // void sendTileData(uint64_t cid, std::vector<int> dests) {
  //   auto c = getTileData(cid);

  //   size_t j = sent_tile_messages.size();

  //   for (auto dest: dests) {
  //     MPI_Request req;
  //     sent_tile_messages.push_back( req );

  //     // TODO: use internal tile members as
  //     c.send(dest);
  //     // c has a switch variable to select MPI msg state
  //
  //     MPI_Isend(
  //         c,
  //         1,
  //         mpi_tile_t,
  //         dest,
  //         commType::CELLDATA,
  //         comm,
  //         &sent_tile_messages[j]
  //         );
  //     j++;
  //   }
  // }


  // /// Receive incoming stuff
  // void communicateRecvTiles() {

  //   recv_info_messages.clear();
  //   recv_tile_messages.clear();

  //   size_t i = 0;
  //   for (int source=0; source<Nrank; source++) {
  //     if (source == rank) { continue; } // do not receive from myself

  //     // communicate with how many tiles there are incoming

  //     // TODO: use MPI_IProbe to check if there are 
  //     // any messages for me instead of assuming that there is

  //     MPI_Request req;
  //     recv_info_messages.push_back( req );

  //     uint64_t Nincoming_tiles;
  //     MPI_Irecv(
  //         &Nincoming_tiles,
  //         1,
  //         MPI_UNSIGNED_LONG_LONG,
  //         source,
  //         commType::NCELLS,
  //         comm,
  //         &recv_info_messages[i]
  //         );

  //     // TODO: Remove this code block and do in background instead
  //     MPI_Wait(&recv_info_messages[i], MPI_STATUS_IGNORE);

  //     /*
  //        fmt::print("{}: I got a message! Waiting {} tiles from {}\n",
  //        rank, Nincoming_tiles, source);
  //        */


  //     // Now receive the tiles themselves
  //     size_t j = recv_tile_messages.size();
  //     for (size_t ic=0; ic<Nincoming_tiles; ic++) {
  //       Tile inc_c(0,0,0, Nx, Ny); // TODO: initialize with better default values

  //       MPI_Request reqc;
  //       recv_tile_messages.push_back( reqc );
  //       MPI_Irecv(
  //           &inc_c,
  //           1,
  //           mpi_tile_t,
  //           source,
  //           commType::CELLDATA,
  //           comm,
  //           &recv_tile_messages[j]
  //           );

  //       MPI_Wait(&recv_tile_messages[j], MPI_STATUS_IGNORE);
  //       j++;

  //       uint64_t cid = inc_c.cid;
  //       if (this->tiles.count(cid) == 0) {
  //         // Tile does not exist yet; create it
  //         // TODO: Check validity of the tile better
  //         // TODO: Use add_tile() instead of directly 
  //         //       probing the class interiors

  //         inc_c.local = false;
  //         tiles.insert( std::make_pair(cid, &inc_c) );

  //       } else {
  //         // Tile is already on my virtual list; update
  //         // TODO: use = operator instead.
  //         auto c = getTileData(cid);

  //         if (c->local) {
  //           // TODO: better error handling; i.e. resolve the conflict
  //           // TODO: use throw exceptions
  //           // fmt::print("{}: ERROR trying to add virtual tile that is already local\n", rank);
  //           exit(1);
  //         }

  //         c->owner             = inc_c.owner;
  //         c->i                 = inc_c.i;
  //         c->j                 = inc_c.j;
  //         c->top_virtual_owner = inc_c.top_virtual_owner;
  //         c->communications    = inc_c.communications;
  //         c->number_of_virtual_neighbors = inc_c.number_of_virtual_neighbors;


  //       };

  //     }
  //     i++;
  //   }
  // }


}; // end of Node class

} // end of corgi namespace


