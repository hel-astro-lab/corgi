#pragma once

// #include <fmt/format.h>
// #include <fmt/format.cc>
// #include <fmt/string.h>
// #include <fmt/ostream.h>

#include <vector>
#include <set>
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
  corgi::tools::sparse_grid<int, D> _mpi_grid;

  /// global large scale block grid where load balance 
  //information is stored
  corgi::tools::sparse_grid<double, D> _work_grid;

  // --------------------------------------------------
  private:
    
  // Mappings
  using TileID_t = uint64_t;
  using Tile_t   = corgi::Tile<D>;
  using Tileptr  = std::shared_ptr<Tile_t>;
  using Tile_map = std::unordered_map<TileID_t, Tileptr>;



  public:

  /// Map with tile_id & tile data
  Tile_map tiles;


  public:
  // --------------------------------------------------
  // Python bindings for mpi_grid

  // get element
  template<typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value, int > 
  py_get_mpi_grid(Indices... indices)  /*const*/
  {
    return _mpi_grid(indices...);
  }


  // set element
  template<typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value, void > 
  py_set_mpi_grid(int val, Indices... indices) {
    _mpi_grid(indices...) = val;
  }


  // get element
  template<typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value, double > 
  py_get_work_grid(Indices... indices)  /*const*/
  {
    return _work_grid(indices...);
  }

  // set element
  template<typename... Indices>
  corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
  corgi::internals::are_integral<Indices...>::value, void > 
  py_set_work_grid(double val, Indices... indices) {
    _work_grid(indices...) = val;
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
  Node(DimensionLength... dimension_lengths) :
    _lengths {{static_cast<size_type>(dimension_lengths)...}},
    _mpi_grid(dimension_lengths...),
    _work_grid(dimension_lengths...),
    env(),
    comm()
  { }


  /*
   * try specializing handy shortcuts to symmetrize construction always assuming 3D input
  template< typename = corgi::internals::enable_if_t< (D == 1), void > > 
  Node(size_t i, size_t j, size_t k) :
    _lengths {{i}},
    _mpi_grid({{i}})
  { }

  template< typename = corgi::internals::enable_if_t< (D == 2), void > > 
  Node(size_t i, size_t j, size_t k) :
    _lengths {{i, j}},
    _mpi_grid({{i, j}})
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
   *  coeffs[i] = ct_accumulate(dimension_lengths, i + 1, Dimensions - i - 1,
   *                                        static_cast<size_type>(1),
   *                                        ct_prod<size_type>);
   *
   */
  std::array<size_type, D>
  compute_index_coeffs(const ::std::array<size_type, D>& dimension_lengths) const noexcept
  {
      std::array<size_type, D> coeffs;
      for (size_type i = 0; i < D; ++i)
      {
          coeffs[i] = corgi::internals::ct_accumulate(
              dimension_lengths,
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

  template <size_t... Is>
  index_type id_impl(
      const corgi::internals::tuple_of<D, size_t>& tuple, 
      std::index_sequence<Is...>) const
  {
    return id( std::get<Is>(tuple)... );
  }

  /// unpack tuple into variadic argument list
  template<typename Indices = std::make_index_sequence<D>>
  index_type id( corgi::internals::tuple_of<D, size_t>& indices)
  {
      return id_impl(indices, Indices{} );
  }

  template<typename Indices = std::make_index_sequence<D>>
  index_type id( const corgi::internals::tuple_of<D, size_t>& indices) const
  {
      return id_impl(indices, Indices{} );
  }
  
  /// Inverse Morton Z-code
  //
  // TODO: make N-dimensional
  corgi::internals::tuple_of<1, index_type> id2index(
      uint64_t cid, 
      std::array<size_type,1> /*lengths*/)
  {
    corgi::internals::tuple_of<1, index_type> indices = std::make_tuple(cid);

    return indices;
  }

  corgi::internals::tuple_of<2, index_type> id2index(
      uint64_t cid,
      std::array<size_type,2> lengths)
  {
    corgi::internals::tuple_of<2, index_type> indices = std::make_tuple
      (
       cid % lengths[0],
      (cid / lengths[0]) % (lengths[1] )
       );

    return indices;
  }

  corgi::internals::tuple_of<3, index_type> id2index(
      uint64_t cid,
      std::array<size_type,3> lengths)
  {
    corgi::internals::tuple_of<3, index_type> indices = std::make_tuple
      (
       cid % lengths[0],
      (cid / lengths[0]) % (lengths[1] ),
       cid /(lengths[0] * lengths[1] )
       );

    return indices;
  }


  public:

  // --------------------------------------------------
  // apply SFINAE to create some shortcuts (when appropriate) 
  // NOTE: valid up to D=3 with x/y/z

  // return global grid sizes
  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=1), T> 
  get_Nx() { return _lengths[0]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=2), T> 
  get_Ny() { return _lengths[1]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=3), T> 
  get_Nz() { return _lengths[2]; }


  // return global grid limits
  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=1), T> 
  get_xmin() { return _mins[0]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=2), T> 
  get_ymin() { return _mins[1]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=3), T> 
  get_zmin() { return _mins[2]; }


  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=1), T> 
  get_xmax() { return _maxs[0]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=2), T> 
  get_ymax() { return _maxs[1]; }

  template<typename T = size_type>
  corgi::internals::enable_if_t< (D>=3), T> 
  get_zmax() { return _maxs[2]; }


  /// Set physical grid size
  void set_grid_lims(
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
  // void add_tile(Tile& tile) {
  //
  // FIXME
  void add_tile(
    Tileptr tileptr,
    corgi::internals::tuple_of<D, size_t> indices
    )
  {

    // claim unique ownership of the tile (for unique_ptr)
    // std::unique_ptr<corgi::Tile> tileptr = std::make_unique<corgi:Tile>(tile);
    // Tileptr tileptr = std::make_unique<Tile_t>(tile);
    
    // calculate unique global tile ID
    uint64_t cid = id( indices );

    // Erase any existing tiles to avoid emplace of doing nothing TODO: is this correct?
    tiles.erase(cid);

    tileptr->index               = indices;
    tileptr->cid                 = cid;
    tileptr->communication.cid   = cid;
    tileptr->communication.owner = comm.rank();
    //tileptr->communication.local = true; //TODO Catch error if tile is not already mine?
    tileptr->lengths             = _lengths;

    // copy indices from tuple into D=3 array in Communication obj
    auto tmp = corgi::internals::into_array(indices);
    for(size_t i=0; i<D; i++) tileptr->communication.indices[i] = tmp[i];

    // tiles.emplace(cid, std::move(tileptr)); // unique_ptr needs to be moved
    tiles.emplace(cid, tileptr); // NOTE using c++14 emplace to avoid copying
    //tiles.insert( std::make_pair(cid, tileptr) ); // NOTE using c++14 emplace to avoid copying
      
    // add to my internal listing
    _mpi_grid( indices ) = comm.rank();
  }


  /// Replace/add incoming tile
  //
  // FIXME
  void replace_tile(
    Tileptr tileptr,
    corgi::internals::tuple_of<D, size_t> indices
    )
  {
    // calculate unique global tile ID
    uint64_t cid = id(indices);

    // add tile if it does not exist
    if(tiles.count(cid) == 0) return add_tile(tileptr, indices);

    // else replace previous one; copy Communication object
    auto& tile = get_tile(cid);
    auto cm = tile.communication;

    tileptr->index   = indices;
    tileptr->cid     = cid;
    tileptr->lengths = _lengths;

    tiles.erase(cid);
    tiles.emplace(cid, tileptr); 

    update_tile(cm);
  }


  /// Shortcut for creating raw tiles with only the internal meta info.
  // to be used with message passing (w.r.t. add_tile that is for use with initialization)
  // FIXME
  void create_tile(Communication& cm)
  {
    //m_author = std::make_shared<Author>(t_author);
    auto tileptr = std::make_shared<Tile_t>();
    tileptr->load_metainfo(cm);

    // additional node info
    tileptr->lengths = _lengths;
    // owner
    // local

    // add
    tiles.emplace(cm.cid, tileptr); // NOTE using c++14 emplace to avoid copying
    _mpi_grid( tileptr->index ) = cm.owner;
  }

  /// Update tile metadata
  // FIXME
  void update_tile(Communication& cm)
  {
    auto& tile = get_tile(cm.cid);
    tile.load_metainfo(cm);
    _mpi_grid( tile.index ) = cm.owner;
  }


  /*! Return a vector of tile indices that fulfill a given criteria.  */
  std::vector<uint64_t> get_tile_ids(
      const bool sorted=false ) 
  {
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
   * Other_t& get_othertype(const std::string& name)
   * {
   *     auto it = otMap.find(name);
   *     if (it == otMap.end()) throw std::invalid_argument("entry not found");
   *     return *(it->second);
   * }
   *
   * This way map retains its ownership of the tile and we avoid giving pointers
   * away from the Class.
   */
  Tile_t& get_tile(const uint64_t cid) {
    auto it = tiles.find(cid);
    if (it == tiles.end()) { throw std::invalid_argument("tile entry not found"); }

    return *(it->second);
  }

  template<typename... Indices>
    corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
    corgi::internals::are_integral<Indices...>::value, 
  Tile_t&>
  get_tile_ind(const Indices... indices) {
    uint64_t cid = id(indices...);
    return get_tile(cid);
  }

  /// \brief Get individual tile (as a pointer)
  Tileptr get_tileptr(const uint64_t cid) {
    auto it = tiles.find(cid);
    //if (it == tiles.end()) { throw std::invalid_argument("entry not found"); }
    if (it == tiles.end()) { return nullptr; };
    return it->second;
  }

  template<typename... Indices>
    corgi::internals::enable_if_t< (sizeof...(Indices) == D) && 
    corgi::internals::are_integral<Indices...>::value, 
  Tileptr>
  get_tileptr_ind(const Indices... indices)
  {
    uint64_t cid = id(indices...);
    return get_tileptr(cid);
  }

  Tileptr get_tileptr(const std::tuple<size_t> ind) {
    size_t i = std::get<0>(ind);
    return get_tileptr_ind(i);
  }

  Tileptr get_tileptr(const std::tuple<size_t, size_t> ind) {
    size_t i = std::get<0>(ind);
    size_t j = std::get<1>(ind);
    return get_tileptr_ind(i, j);
  }

  Tileptr get_tileptr(const std::tuple<size_t, size_t, size_t> ind) {
    size_t i = std::get<0>(ind);
    size_t j = std::get<1>(ind);
    size_t k = std::get<2>(ind);
    return get_tileptr_ind(i, j, k);
  }


  /// Return all local tiles
  std::vector<uint64_t> get_local_tiles(
      const bool sorted=false ) {

    std::vector<uint64_t> tile_list = get_tile_ids(sorted);
    std::vector<uint64_t> ret;
    ret.reserve(tile_list.size());

    for(auto elem : tile_list) {
      if(tiles.at( elem )->communication.owner == comm.rank() ) {
        ret.push_back(elem);
      }
    }
    return ret;

    //size_t i = 0, len = tile_list.size();
    //while(i < len) {
    //  std::cout << comm.rank() << ": " << i << "/" << len << "\n";
    //  if(!tiles.at( tile_list[i] )->communication.local) {
    //    std::swap(tile_list[i], tile_list.back());
    //    tile_list.pop_back();
    //    len--;
    //  } else {
    //    i++;
    //  }
    //}
    //return tile_list;
  }


  /// Return all tiles that are of VIRTUAL type.
  std::vector<uint64_t> get_virtuals(
      const bool sorted=false ) 
  {

    std::vector<uint64_t> tile_list = get_tile_ids(sorted);
    std::vector<uint64_t> ret;
    ret.reserve(tile_list.size());

    for(auto elem : tile_list) {
      if(tiles.at( elem )->communication.owner != comm.rank() ) {
        ret.push_back(elem);
      }
    }
    return ret;

    //size_t i = 0, len = tile_list.size();
    //while(i < len) {
    //  if (tiles.at( tile_list[i] )->communication.local) {
    //    std::swap(tile_list[i], tile_list.back());
    //    tile_list.pop_back();
    //    len--;
    //  } else {
    //    i++;
    //  }
    //}
    //return tile_list;
  }

  /// Return all local boundary tiles
  // TODO: update this to use the internal boundary_tile_map
  std::vector<uint64_t> get_boundary_tiles(
      const bool sorted=false ) {

    std::vector<uint64_t> tile_list = get_tile_ids(sorted);

    size_t i = 0, len = tile_list.size();
    while(i < len) {

      // remove if there are no virtual nbors and tile is not mine -> opposite means its boundary
      if (tiles.at( tile_list[i] )-> communication.number_of_virtual_neighbors == 0 || 
          tiles.at( tile_list[i] )-> communication.owner != comm.rank()
          ) {
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
  bool is_local(uint64_t cid) {
    bool local = false;

    // Do we have it on the list?
    if (tiles.count( cid ) > 0) {
      // is it local (i.e., not virtual)
      if ( tiles.at(cid)->communication.owner == comm.rank() ) {
        local = true;
      }
    }

    return local;
  }

  /// return all virtual tiles around the given tile
  std::vector<uint64_t> virtual_nhood(uint64_t cid) 
  {
    auto& c = get_tile(cid);
    auto neigs = c.nhood();

    std::vector<uint64_t> vnhood;
    for(auto& indx: neigs) {
      if(_mpi_grid(indx) != comm.rank()) vnhood.push_back( id(indx) );
    }

    return vnhood;
  }

  /// return all owners of virtual tiles around the given tile
  std::vector<int> virtual_nhood_owners(uint64_t cid) 
  {
    auto& c = get_tile(cid);
    auto neigs = c.nhood();

    std::vector<int> virtual_owners;
    for(auto& indx: neigs) {
      int whoami = _mpi_grid(indx); // Get tile id from index notation
      if(whoami != comm.rank()) virtual_owners.push_back( whoami );
    }

    return virtual_owners;
  }


  /// map of owners to exterior (=virtual) tiles
  std::map<int, std::set<uint64_t> > virtual_tile_list;

  /// map of my (local) tiles to exterior ranks
  std::map<uint64_t, std::set<int> > boundary_tile_list;


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
  void analyze_boundaries() {

    virtual_tile_list.clear();
    boundary_tile_list.clear();

    // analyze all of my local tiles
    for(auto cid: get_local_tiles()) {
      auto& c = get_tile(cid);

      // analyze c's neighborhood
      auto neigs = c.nhood();
      for(auto& indx: neigs) {
        int whoami = _mpi_grid(indx); // Get tile id from index notation

        // if nbor tile is virtual
        if(whoami != comm.rank()) {

          // then the local tile cid is a boundary tile to this rank (=whoami)
          boundary_tile_list[cid].insert(whoami);

          // and then ncid is my virtual exterior tile
          uint64_t ncid = id(indx);
          virtual_tile_list[whoami].insert(ncid);
        }
      }

      // mark completely local if no virtuals around this tile
      // overwritten in next loop based on real values
      // could also carry a flag through the loop above to check this...
      c.communication.number_of_virtual_neighbors = 0;
      c.communication.communications              = 0;
    }

    // update this info into my local boundary tiles
    for(auto&& elem : boundary_tile_list) {
      auto& c = get_tile(elem.first);

      // set into vector
      std::vector<int> virtual_owners(elem.second.begin(), elem.second.end()); 

      //c.communication.top_owner = top_owner; // can not be computed from set
      c.communication.communications = virtual_owners.size();
      c.communication.number_of_virtual_neighbors = virtual_owners.size(); // not same as original
      c.communication.virtual_owners = virtual_owners;

      // add to send queue
      uint64_t cid = elem.first;
      if (std::find( send_queue.begin(), send_queue.end(), cid) 
            == send_queue.end()) {
          send_queue.push_back( cid );
          send_queue_address.push_back( virtual_owners );
      }
    }

    //TODO: can also pre-create virtual tiles (if not existing in node yet)  

  }


  //void analyze_boundaries_old() {

  //  // old tile information update  
  //  for(auto cid: get_local_tiles()) {
  //    //std::cout << comm.rank() << ": ab: " << cid << "\n";
  //    std::vector<int> virtual_owners = virtual_nhood_owners(cid);
  //    size_t N = virtual_owners.size();
  //    //std::cout << comm.rank() << ": ab: " << cid << "N:" << N << "\n";

  //    // If N > 0 then this is a local boundary tile.
  //    // other criteria could also apply but here we assume
  //    // neighborhood according to spatial distance.
  //    if (N > 0) {

  //      /* Now we analyze `owner` vector as:
  //       * - sort the vector
  //       * - compute mode of the list to see who owns most of them
  //       * - remove repeating elements creating a unique list. */

  //      // sort
  //      std::sort( virtual_owners.begin(), virtual_owners.end() );

  //      // compute mode by creating a frequency array
  //      // NOTE: in case of same frequency we implicitly pick smaller rank
  //      int max=0, top_owner = virtual_owners[0];
  //      for(size_t i=0; i<virtual_owners.size(); i++) {
  //        int co = (int)count(virtual_owners.begin(), 
  //            virtual_owners.end(), 
  //            virtual_owners[i]);
  //        if(co > max) {      
  //          max = co;
  //          top_owner = virtual_owners[i];
  //        }
  //      } 

  //      // remove duplicates
  //      virtual_owners.erase( unique( virtual_owners.begin(), 
  //            virtual_owners.end() 
  //            ), virtual_owners.end() );


  //      // update tile values
  //      auto& c = get_tile(cid);
  //      c.communication.top_virtual_owner = top_owner;
  //      c.communication.communications    = virtual_owners.size();
  //      c.communication.number_of_virtual_neighbors = N;
  //      c.communication.virtual_owners    = virtual_owners;

  //      if (std::find( send_queue.begin(), send_queue.end(), cid) 
  //          == send_queue.end()) 
  //      {
  //        send_queue.push_back( cid );
  //        send_queue_address.push_back( virtual_owners );
  //      }



  //    } else { // else N == 0
  //      auto& c = get_tile(cid);
  //      c.communication.number_of_virtual_neighbors = 0;
  //      c.communication.communications              = 0;
  //    }

  //  }
  //}


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

  std::vector<mpi::request> sent_info_messages;
  std::vector<mpi::request> sent_tile_messages;
  std::unordered_map<int, std::vector<mpi::request>> sent_data_messages;
  //std::unordered_map<int, std::vector<MPI_Request>> sent_data_messages;

  std::vector<mpi::request> recv_info_messages;
  std::vector<mpi::request> recv_tile_messages;
  std::unordered_map<int, std::vector<mpi::request>> recv_data_messages;
  //std::unordered_map<int, std::vector<MPI_Request>> recv_data_messages;

  std::vector<mpi::request> sent_adoption_messages;
  std::vector<mpi::request> recv_adoption_messages;

  // /// Broadcast master ranks mpi_grid to everybody
  void bcast_mpi_grid() {

    // total size
    int N = 1;
    for (size_t i = 0; i<D; i++) N *= _lengths[i];
    std::vector<int> tmp;

    if (comm.rank() == 0) {
      tmp = _mpi_grid.serialize();
    } else {
      tmp.resize(N);
      for(int k=0; k<N; k++) {tmp[k] = -1;};
    }

    MPI_Bcast(&tmp[0],
        N, 
        MPI_INT, 
        0, 
        MPI_COMM_WORLD
        );

    // unpack
    if(comm.rank() != 0) {
      _mpi_grid.deserialize(tmp, _lengths);
    }
  }

  /// update work arrays from other nodes and send mine
  void allgather_work_grid() 
  {

    // total size
    int N = 1;
    for (size_t i = 0; i<D; i++) N *= _lengths[i];

    // message buffers
    std::vector<double> recv(N*comm.size());
    std::vector<double> orig = _work_grid.serialize();

    // mask all work values that are not mine
    std::vector<int> ranks   = _mpi_grid.serialize();
    for(size_t i=0; i<ranks.size(); i++) {
      if( ranks[i] != comm.rank() ) orig[i] = -1.0;
    }

    MPI_Allgather(
        &orig[0],
        orig.size(), 
        MPI_DOUBLE, 
        &recv[0], 
        orig.size(),
        MPI_DOUBLE,
        MPI_COMM_WORLD
        );
    
  
    std::vector<double> new_work(N);
    std::fill(new_work.begin(), new_work.end(), -1.0);

    double val;
    for(int j=0; j<comm.size(); j++) {
      for(int i=0; i<N; i++) {
        val = recv[j*N + i];
        if(val != -1.0) {
          assert(new_work[i] == -1.0); // test that we dont overwrite
          new_work[i] = val;
        }
      }
    }
      
    // test that we did not miss elements
    for(auto& val : new_work) assert(val != -1.0); 

    // upload back to grid
    _work_grid.deserialize(new_work, _lengths);

  }


  // Update work load grid from my local tiles
  void update_work()
  {
    for(auto& cid : get_local_tiles()) {
      auto& tile = get_tile(cid);
     _work_grid( tile.index ) = tile.get_work();
    }
  }



  /// Issue isends to everywhere
  // First we send a warning message of how many tiles to expect.
  // Based on this the receiving side can prepare accordingly.
  void send_tiles() {

    sent_info_messages.clear();
    sent_tile_messages.clear();

    //for (int dest = 0; dest<comm.size(); dest++) {
    //  if( dest == comm.rank() ) { continue; } // do not send to myself

    //  int i = 0;
    //  std::vector<int> to_be_sent;
    //  for(std::vector<int> address: send_queue_address) {
    //    if( std::find( address.begin(),
    //          address.end(),
    //          dest) != address.end()) 
    //    {
    //      to_be_sent.push_back( i );
    //    }
    //    i++;
    //  }

    //  // initial message informing how many tiles are coming
    //  // TODO: this whole thing could be avoided by using 
    //  // MPI_Iprobe in the receiving end. Maybe...
    //  auto number_of_incoming_tiles = static_cast<int>(to_be_sent.size());

    //  //std::cout << comm.rank() 
    //  //          << " sending message to " 
    //  //          << dest
    //  //          << " incoming number of tiles " 
    //  //          << number_of_incoming_tiles
    //  //          << "\n";

    //  mpi::request req;
    //  req = comm.isend(dest, commType::NTILES, number_of_incoming_tiles);
    //  sent_info_messages.push_back( req );

    //}

    // send the real tile meta info data now
    // We optimize this by only packing the tile data
    // once, and then sending the same thing to everybody who needs it.
    // FIXME: not really...
    //int i = 0;
    //for(auto cid: send_queue) {
    //  auto& tile = get_tile(cid);
    //  for(int dest: send_queue_address[i]) {

    //    mpi::request req;
    //    req = comm.isend(dest, commType::TILEDATA, tile.communication);

    //    sent_tile_messages.push_back( req );
    //  }
    //  i++;
    //}

    // send all tiles
    for(auto&& elem : boundary_tile_list) {
      auto& tile = get_tile(elem.first);

      std::cout << comm.rank() << " sending cid message " << elem.first << " ---> ";

      for(int dest : elem.second) {
        std::cout << "," << comm.rank() << ":" << dest;

        mpi::request req;
        req = comm.isend(dest, commType::TILEDATA, tile.communication);

        sent_tile_messages.push_back( req );
      }

      std::cout << "\n";
    }

  }

  /// Send individual tile to dest
  void send_tile(uint64_t cid, int dest)
  {
    mpi::request req;

    auto& tile = get_tile(cid);
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
    //rcom.local = false; // received tiles are automatically virtuals
    create_tile(rcom);

    return;
  }

  /// Receive incoming stuff
  //void recv_tiles() {

  //  recv_info_messages.clear();
  //  recv_tile_messages.clear();

  //  //size_t i = 0;
  //  for (int source=0; source<comm.size(); source++) {
  //    if (source == comm.rank() ) continue; // do not receive from myself

  //    // communicate with how many tiles there are incoming

  //    // TODO: use MPI_IProbe to check if there are 
  //    // any messages for me instead of assuming that there is

  //    // TODO: encapsulate into vector that can be received & 
  //    // processed more later on

  //    int number_of_incoming_tiles=0;
  //    mpi::request req;
  //    req = comm.irecv(source, commType::NTILES, number_of_incoming_tiles);

  //    // TODO: Remove this code block and do in background instead
  //    req.wait();
  //    //recv_info_messages.push_back( req );

  //    /*
  //       fmt::print("{}: I got a message! Waiting {} tiles from {}\n",
  //       rank, number_of_incoming_tiles, source);
  //       */
  //    std::cout << comm.rank()
  //              << " I got a message! Waiting " 
  //              << number_of_incoming_tiles << " tiles from " 
  //              << source
  //              << "\n";

  //    // Now receive the tiles themselves
  //    for (int ic=0; ic<number_of_incoming_tiles; ic++) {
  //      mpi::request reqc;
  //      Communication rcom;

  //      reqc = comm.irecv(source, commType::TILEDATA, rcom);
  //      
  //      // TODO non blocking
  //      reqc.wait();
  //      recv_tile_messages.push_back( reqc );
  //        
  //      if(this->tiles.count(rcom.cid) == 0) { // Tile does not exist yet; create it
  //        // TODO: Check validity of the tile better
  //        create_tile(rcom);
  //      } else { // Tile is already on my virtual list; update
  //        update_tile(rcom);  
  //      };

  //      // update and make non-local virtual
  //      //auto& new_tile = get_tile(rcom.cid);
  //      //new_tile.communication.local = false;
  //    }
  //  }

  //  // process all mpi requests; otherwise we leak memory
  //  mpi::wait_all(recv_info_messages.begin(), recv_info_messages.end());
  //  mpi::wait_all(sent_info_messages.begin(), sent_info_messages.end());

  //  mpi::wait_all(recv_tile_messages.begin(), recv_tile_messages.end());
  //  mpi::wait_all(sent_tile_messages.begin(), sent_tile_messages.end());



  //  clear_send_queue();
  //}


  std::vector<Communication> rcoms;
  void recv_tiles() {
    recv_tile_messages.clear();

    int i = 0;
    for(auto&& elem : virtual_tile_list) {
      int orig = elem.first;

      std::cout << comm.rank() << " receiving from "<<  elem.first << " <--- ";
      for(uint64_t cid : elem.second) {
        mpi::request reqc;

        //rcoms.emplace_back();
        Communication rcom;
        rcoms.push_back(rcom);
        std::cout << "," << comm.rank() << ":" << cid << "(" << rcoms.size()<<")" ;
        reqc = comm.irecv(orig, commType::TILEDATA, rcoms.at(i) );
        
        recv_tile_messages.push_back( reqc );
        i++;
      }

      std::cout << "\n";
    }

    std::cout << comm.rank() << " waiting...\n";

    // process all mpi requests; otherwise we leak memory
    mpi::wait_all(recv_tile_messages.begin(), recv_tile_messages.end());
    std::cout << comm.rank() << " unpacking...\n";

    // unpack here
    for(auto rcom : rcoms) {
      if(this->tiles.count(rcom.cid) == 0) { // Tile does not exist yet; create it
        // TODO: Check validity of the tile better
        create_tile(rcom);
      } else { // Tile is already on my virtual list; update
        update_tile(rcom);  
      }
    }

    std::cout << comm.rank() << " waiting sents...\n";

    // wait rest of messages too
    mpi::wait_all(sent_tile_messages.begin(), sent_tile_messages.end());

    std::cout << comm.rank() << " done with sents...\n";
  }


  // Decide who to adopt
  //
  // Every node has their own council, but they should all come to same 
  // conclusion.
  
  private:
  std::vector<int> adoptions; 
  std::vector<int> kidnaps; 
  double min_quota =  0.05;
  double max_quota =  0.05; // in fraction of all tiles per color



  public:
  /// Compute maximum number of new tiles I can adopt
  double get_quota(int rank)
  {
    // integrate over all work
    int N = 1;
    for (size_t i = 0; i<D; i++) N *= _lengths[i];
    
    // ideal work balance
    double ideal_work = 0.0;
    for(auto& elem : _work_grid) ideal_work += elem.second;
    ideal_work /= comm.size();

    // current work load
    double current_workload = 0.0;
    for(auto& elem : _mpi_grid) {
      if(elem.second == rank) {
        current_workload += _work_grid( elem.first );
      }
    }

    /// excess work I can do
    double excess = ideal_work - current_workload;
    //excess = excess > 0.0 ? excess : 0.0;


    //if(rank == comm.rank()) {
    //std::cout << comm.rank() << " get_quota for rank gives "
    //  << "ideal_work:" << ideal_work
    //  << "current   :" << current_workload
    //  << "excess   :" << excess
    //  << "\n";
    //}

    return excess;

    double quota = excess > max_quota*ideal_work ? max_quota*ideal_work : excess;
    return quota > -min_quota*ideal_work ? quota : -min_quota*ideal_work;
  }
  
  /// Propagate CA rules one step forward and decide who adopts who
  void adoption_council()
  {
    int quota = (int)get_quota(comm.rank());

    // collect virtual tiles and their metainfo into a container
    std::vector<Communication> virtuals;
    for(auto cid : get_virtuals() ) {
      auto& tile = get_tile(cid);
      virtuals.push_back(tile.communication);
    }

    // sort in descending order using lambda function
    std::sort(virtuals.begin(), virtuals.end(), 
        [](const auto& lhs, const auto& rhs) 
    {
      return lhs.number_of_virtual_neighbors > rhs.number_of_virtual_neighbors;
    });


    // now loop over all virtual tiles and check if I can adopt someone
    adoptions.clear();
    for(auto& vir : virtuals) {

      //if(vir.number_of_virtual_neighbors <= 3) continue;

      // skip tiles where I am not the top owner
      if(vir.top_virtual_owner == comm.rank()) {
        adoptions.push_back( vir.cid );

        //std::cout << comm.rank() 
        //  << ": adoption council marks " << vir.cid
        //  << " at " << vir.indices[0] << " " << vir.indices[1] << "\n";
      }

      if( (int)adoptions.size() >= quota) break;
    }

  }

  /// general N-dim implementation of wrap
  size_t wrap(int ind, size_t d)
  {
    auto N = static_cast<int>(_lengths[d]);
    while (ind < 0) { ind += N; }
    while (ind >= N) { ind -= N; }
    return static_cast<size_t>(ind);
  }


  /// return index of tiles in relative to my position
  const corgi::internals::tuple_of<D, size_t> neighs(
      corgi::internals::tuple_of<D, int> indices,
      corgi::internals::tuple_of<D, int> indices_rel)
  {
    auto cur = corgi::internals::into_array(indices);
    auto rel = corgi::internals::into_array(indices_rel);

    for(size_t i=0; i<D; i++) {
      cur[i] = static_cast<size_t>(
          wrap( static_cast<int>(rel[i]) + 
            static_cast<int>(cur[i]), i)
          );
    }

    auto ret = corgi::internals::into_tuple(cur);
    return ret;
  }

  /// Return full Moore neighborhood around me
  std::vector< corgi::internals::tuple_of<D, size_t> > nhood(
      corgi::internals::tuple_of<D, size_t> indices)
  {
    // FIXME; generalize the size (now assumes D=2)
    std::vector< corgi::internals::tuple_of<D, size_t> > nh;
    nh.reserve(8);

    for(auto& reli : corgi::ca::moore_neighborhood<D>() ){
      nh.push_back( neighs(indices, reli) );
    }
    return nh;
  }


  void adoption_council2()
  {
    adoptions.clear();
    kidnaps.clear();
    std::vector<double> alives(comm.size());

    //int myrank = comm.rank();

    // for updated values
    corgi::tools::sparse_grid<int, D> new_mpi_grid(_mpi_grid);

    // radius of Gaussian kernel
    int dt = sqrt(*std::max_element(_lengths.begin(), _lengths.end() )); 
    //int dt = *std::max_element(_lengths.begin(), _lengths.end() )/2; 
    //dt += 0.1*dt;

    // gaussian kernel; i.e., relative indices how we convolve
    auto kernel = corgi::ca::chessboard_neighborhood<D>(dt);


    // keep track of tile changes
    std::vector<int> obtained(comm.size()), lost(comm.size());
    std::fill(obtained.begin(), obtained.end(), 0); // reset vector
    std::fill(lost.begin(), lost.end(), 0); // reset vector

    // transferred work
    std::vector<double> transfers(comm.size());
    std::fill(transfers.begin(), transfers.end(), 0.0); // reset vector

    // absolute quota in units of work
    std::vector<double> quota(comm.size());
    for(int rank=0; rank<comm.size(); rank++) quota[rank] = get_quota(rank);

    // relative quota
    std::vector<double> rel_quota(comm.size());
    double total_work = 0.0;
    for(auto& elem : _work_grid) total_work += elem.second;
    for(size_t i=0; i<rel_quota.size(); i++) rel_quota[i] = comm.size()*quota[i]/total_work;

    //std::cout << comm.rank() << ": my quota : " << quota[myrank] 
    //  << " (" << rel_quota[myrank] << ")"
    //  << "\n";

    //--------------------------------------------------

    // process the complete grid (including remote neighbors)
    int new_color;
    for(auto&& elem : _mpi_grid) {
      auto& ind     = elem.first;
      int old_color = elem.second;

      std::fill(alives.begin(), alives.end(), 0.0); // reset vector

      // resolve neighborhood; diffusion step
      //alives[old_color] = 1.0;
      alives[old_color] = (1.0/4.0/M_PI/static_cast<double>(dt));
      //alives[old_color] = sqrt(0.5/M_PI);

      // Limited Moore nearest neighborhood
      //auto neigs = nhood(ind);
      //for(auto& nindx : neigs) {
      //  color = _mpi_grid(nindx);
      //  assert(0 <= color && color < comm.size());

      //  alives[color]++; 
      //}


      // full Gaussian kernel
      double r;
      int color;
      for(auto& reli : kernel ) {
        auto nindx = neighs(ind, reli);
        color = _mpi_grid(nindx);

        //r = geom::eulerian_distance(reli);  
        r = geom::manhattan_distance<D>(reli);  
        //r = geom::chessboard_distance<D>(reli);  

        //alives[color] += exp(-r*r/static_cast<double>(dt));
        //alives[color] += r/(2.0*M_PI*static_cast<double>(dt));
          
        alives[color] += (1.0/4.0/M_PI/static_cast<double>(dt))
                         *exp(-r*r/(4.0*static_cast<double>(dt)));
          
        //alives[color] += sqrt(0.5/M_PI)*exp(-r*r/(2.0));

        //alives[color] += exp(-r*r/4.0);
      }


      // normalize
      double norm = 0.0;
      for(auto& val : alives) norm += val;
      for(auto& val : alives) val /= norm;

      // add relative quota for balance 
      for(size_t i=0; i<alives.size(); i++) alives[i] += 0.5*rel_quota[i];


      //auto maxe = std::max_element(alives.begin(), alives.end());
      //if(*maxe >= 0.5/( (double)comm.size() ) ) {
      //  new_color = std::distance(alives.begin(), maxe);
      //} else {
      //  new_color = old_color;
      //}

      // get mode, i.e., most frequent color; sharpening step
      new_color = std::distance( alives.begin(), 
          std::max_element(alives.begin(), alives.end()));

      // stay alive if all values are the same; i.e., we are in equilibrium
      //bool equilibrium = true;
      //int ref_val = alives[0];
      //for(auto& val : alives) {
      //  if(val != ref_val) {
      //    equilibrium = false;
      //    break;
      //  }
      //}
      //if(equilibrium) continue;


      // FIXME
      //if (new_color != old_color) {
      //  uint64_t cid = id(ind);
      //  auto index = id2index(cid, _lengths);
      //  std::cout << comm.rank() << ": tile " << cid 
      //    << " has been kidnapped by evil " << new_color 
      //    << " at (" << std::get<0>(index) << "," 
      //    << " from " << old_color << " --- alives:";
      //  for(auto val : alives) std::cout << " " << val;
      //  double sum = 0.0;
      //  for(auto val : alives) sum += val;
      //  std::cout << " sum: " << sum;
      //  std::cout << "\n";
      //}

      obtained[new_color]++;
      lost[old_color]++;

      // progress one step
      new_mpi_grid(ind) = new_color;
    }

    //std::cout << comm.rank() << ": rank gained/lost= " 
    //  << obtained[myrank] << "/" << lost[myrank] 
    //  << " in += " << obtained[myrank] - lost[myrank] << "\n";


    //--------------------------------------------------


    // next, process changes of ownership
    for(auto&& elem : _mpi_grid) {
      auto& ind     = elem.first;
      int old_color = elem.second;
      int new_color = new_mpi_grid(ind);
      uint64_t cid = id(ind);

      // if no changes, then skip everything
      if(new_color == old_color) continue;


      // check that velocity is not too great; 
      // i.e., change of boundary happens via virtual tiles
      bool is_virtual = false;
      int color;
      for(auto& reli : corgi::ca::moore_neighborhood<D>() ){
        auto nindx = neighs(ind, reli);
        color = _mpi_grid(nindx);
        if(color == new_color) {
          is_virtual = true;
          break;
        }
      }

      // abort if this is not virtual
      if(!is_virtual) {
        //std::cout << comm.rank() << ": XXX tile " << cid 
        //    << " is going too fast to: " << new_color << "\n";
        new_mpi_grid(ind) = old_color;
        continue;
      }
        
      // keep track of work / individual load
      transfers[new_color] += _work_grid(ind);

      // check if we exceed our quota; if yes, then must reinitialize grid
      //if(transfers[new_color] > quota[new_color]) {
      //  new_mpi_grid(ind) = old_color;
      //  transfers[new_color] -= _work_grid(ind);
      //  continue;
      //}

      // I have adopted a tile
      if(new_color == comm.rank()) {
        assert( tiles.count(cid) > 0 );
        adoptions.push_back(cid);

        // FIXME
        //std::cout << comm.rank() 
        //  << ": adoption council marks " << cid << " to be adopted\n";
        //    //<< " at " << tile.communication.indices[0] << " " << tile.communication.indices[1] << "\n";
        auto& tile  = get_tile(cid);
        tile.communication.owner = comm.rank();


        // A tile has been kidnapped from me
      } else if(old_color == comm.rank()) {
        kidnaps.push_back(cid);

        if(is_local(cid)) {
          auto& tile  = get_tile(cid);
          tile.communication.owner = new_color;

          // FIXME
          //std::cout << comm.rank() << ": oh gosh tile " << cid << " is mine! \n";
        } else {
          // something went wrong; I should have this tile but I do not?!
          assert(false);
        }
      }  
      // third option is that somebody else was involved. 

    } // end of loop over elements


    // global progress
    _mpi_grid = std::move(new_mpi_grid);

    // TODO: run analyze_boundaries automatically here to keep 
    //       virtual_tile_list and boundary_tile_list up-to-date.
    //       They are needed e.g. in get_boundary_tiles.

    return;
  }


  /// iterate over tiles in adoption vector and claim them to me
  void adopt()
  {
    for(auto cid : adoptions){
      auto& vir = get_tile(cid);
      //std::cout << comm.rank() 
      //  << ": adopting " << vir.cid
      //  << " at " << vir.communication.indices[0] << " " 
      //            << vir.communication.indices[1] << "\n";


      vir.communication.owner = comm.rank();
      //vir.communication.local = true;

      _mpi_grid( vir.index ) = comm.rank();
    }
  }



  /// send MPI message of my adoptions to everybody
  void send_adoptions()
  {
    sent_adoption_messages.clear();

    // ensure that adoptions vector is of standard length
    //if(adoptions.size() < max_quota) adoptions.resize(max_quota);
    while((int)adoptions.size() < max_quota) adoptions.push_back(-1);

    for (int dest = 0; dest<comm.size(); dest++) {
      if( dest == comm.rank() ) { continue; } // do not send to myself

      mpi::request req;
      req = comm.isend(dest, commType::ADOPT, adoptions.data(), max_quota);
      sent_adoption_messages.push_back( req );
    }
  }

  /// receive MPI adoption messages from others
  void recv_adoptions()
  {
    recv_adoption_messages.clear();

    // ensure that the receiving array is of correct size
    kidnaps.resize(comm.size()*max_quota);

    for (int orig = 0; orig<comm.size(); orig++) {
      if( orig == comm.rank() ) { continue; } // do not recv from myself

      mpi::request req;
      req = comm.irecv(orig, commType::ADOPT, &kidnaps[orig*max_quota], max_quota);
      recv_adoption_messages.push_back( req );
    }

  }

  /// wait and unpack MPI adoption messages
  //
  // Here we try to keep track of all the ownership changes,
  // even if they are remote and do not consider me.
  void wait_adoptions()
  {
    // wait
    mpi::wait_all(recv_adoption_messages.begin(), recv_adoption_messages.end());
    mpi::wait_all(sent_adoption_messages.begin(), sent_adoption_messages.end());

    // unpack
    int kidnapped_cid;
    for (int orig = 0; orig<comm.size(); orig++) {
      if( orig == comm.rank() ) { continue; } // do not process myself

      for(int i=0; i<max_quota; i++) {
        kidnapped_cid = kidnaps[orig*max_quota + i];
        if(kidnapped_cid == -1) continue;

        auto index = id2index(kidnapped_cid, _lengths);

        //std::cout << comm.rank() << ": tile " << kidnapped_cid 
        //  << " has been kidnapped by evil " << orig << "\n";

        if(is_local(kidnapped_cid)) {
          //std::cout << comm.rank() << ": tile " << kidnapped_cid << " is mine! \n";
          //tiles.erase(kidnapped_cid);

          auto& tile = get_tile(kidnapped_cid);
          tile.communication.owner = orig;
          //tile.communication.local = false;
        }

        // update global status irrespective of if it is mine or not
        _mpi_grid(index) = orig;
      }
    }
  }

  /// shortcut for calling blocking version of adoption messaging
  void communicate_adoptions()
  {
    send_adoptions();
    recv_adoptions();
    wait_adoptions();
  }


  /// loop over all virtuals and remove them
  void erase_virtuals()
  {

    int whoami;
    for(auto cid : get_virtuals() ) {

      // check if virtual tile is still needed
      auto& c = get_tile(cid);
      auto neigs = c.nhood();
      for(auto& indx: neigs) {
        whoami = _mpi_grid(indx); 
        if(whoami == comm.rank()) continue;
      }

      tiles.erase(cid);
    }
  }


  // --------------------------------------------------
  // user-data message routines

  /// Initialize vector of vector if needed
  //void initialize_message_array(
  //    std::vector<std::vector<mpi::request>>& arr, 
  //    int tag)
  //{
  //  while((int)arr.size() <= tag) {
  //    std::vector<mpi::request> arri(0);
  //    arr.push_back( arri );
  //  }
  //}

  /// Call mpi send routines from tile for the boundary regions
  // NOTE: we bounce sending back to tile members,
  //       this way methods can be extended for different types of send.
  void send_data(int mode)
  {
    sent_data_messages[mode] = {};

    
    // re-order sends and compute mpi tags
    std::map<int, std::vector<uint64_t> > tags;
    for(auto cid : get_boundary_tiles() ) {
      auto& tile = get_tile(cid);
      for(auto dest: tile.communication.virtual_owners) {
        tags[dest].push_back(cid);
      }
    }
    for(auto& elem : tags) sort(elem.second.begin(), elem.second.end());

    // print
    //for(auto& elem : tags) {
    //  std::cout << comm.rank() << ": S:" << elem.first << " --- ";
    //  for(auto& c : elem.second) std::cout << c << ", ";
    //  std::cout << "\n";
    //}

    int dest;
    uint64_t cid;
    for(auto& elem : tags) {
      dest = elem.first;
      for(int i = 0; i<(int)elem.second.size(); i++) {
        cid = elem.second[i];
        auto& tile = get_tile(cid);
        auto reqs = tile.send_data(comm, dest, mode, i);

        for(auto req : reqs) sent_data_messages.at(mode).push_back(req);
      }
    }

    //for(auto cid : get_boundary_tiles() ) {
    //  auto& tile = get_tile(cid);
    //  for(auto dest: tile.communication.virtual_owners) {
    //    auto reqs = tile.send_data(comm, dest, mode);
    //    
    //    for(auto req : reqs) sent_data_messages.at(mode).push_back(req);
    //  }
    //}
  }


  /// Call mpi recv routines from tile for the virtual regions
  // NOTE: we bounce receiving back to tile members,
  //       this way they can be extended for different types of recv.
  void recv_data(int mode)
  {
    recv_data_messages[mode] = {};

    // re-order sends and compute mpi tags
    std::map<int, std::vector<uint64_t> > tags;
    for(auto cid : get_virtuals() ) {
      auto& tile = get_tile(cid);
      tags[tile.communication.owner].push_back(cid);
    }
    for(auto& elem : tags) sort(elem.second.begin(), elem.second.end());

    // print
    //for(auto& elem : tags) {
    //  std::cout << comm.rank() << ": R:" << elem.first << " --- ";
    //  for(auto& c : elem.second) std::cout << c << ", ";
    //  std::cout << "\n";
    //}


    int orig;
    uint64_t cid;
    for(auto& elem : tags) {
      orig = elem.first;
      for(int i = 0; i<(int)elem.second.size(); i++) {
        cid = elem.second[i];
        auto& tile = get_tile(cid);
        auto reqs = tile.recv_data(comm, orig, mode, i);

        for(auto req : reqs) recv_data_messages.at(mode).push_back(req);
      }
    }

    //int nc = 0, nv = 0;
    //for(auto cid : get_virtuals() ) {
    //  auto& tile = get_tile(cid);
    //  auto reqs = tile.recv_data(comm, tile.communication.owner, mode);

    //  for(auto req : reqs) recv_data_messages.at(mode).push_back(req);
    //  //nc += reqs.size();
    //  //nv++;
    //}

    //std::cout << comm.rank() << ": recv buffer size " << nc << " nv: " << nv << "\n";
  }

  /// barrier until all (primary) data is received
  void wait_data(int tag)
  {
    //assert( tag < (int)recv_data_messages.size() );
    assert( sent_data_messages.count(tag) > 0 );
    assert( recv_data_messages.count(tag) > 0 );
    
    mpi::wait_all( recv_data_messages[tag].begin(), recv_data_messages[tag].end() );
    mpi::wait_all( sent_data_messages[tag].begin(), sent_data_messages[tag].end() );
    //for(auto& req : recv_data_messages[tag]) req.wait();

    // vanilla MPI
    //MPI_Status stats[2];
    //MPI_Waitall( 
    //    recv_data_messages[tag].size(),
    //    &recv_data_messages[tag][0], 
    //    MPI_STATUSES_IGNORE
    //    )
    
    //int n1 = 0;
    //for(auto& vec : sent_data_messages){ n1 += vec.second.size(); }
    //int n2 = 0;
    //for(auto& vec : recv_data_messages){ n2 += vec.second.size(); }
    //std::cout << comm.rank() << ": wait buffer size " << n1 << " " << n2 << "\n";

    // erase (do not force capacity change)
    sent_data_messages[tag] = {};
    recv_data_messages[tag] = {};

    // erase and force clean
    //std::vector<mpi::request>().swap( sent_data_messages[tag] );
    //std::vector<mpi::request>().swap( recv_data_messages[tag] );
    
  }


  void probe_data(int tag)
  {
    // probe for msg
    // put it where?
    // can not get cid
    // where to tmp store it?
    // how to communicate back to tile
    // storage is generally in tile

    return;
  }

  uint64_t reduced_tile_id(uint64_t cid) {
    // alternative is to reduce cid size
    // calc rcid that is used as a tag instead
    // optimal memory usage (no extra transfers)
    // maximum size must be <2^21-1
    // how to calc rcid unambiguously 

    // pre-calc cid & (orig -> dest) number combination?
    // can be done without communications since global knowledge
      
    // algorithm sketch
    // calc neighbors

    //left, right, up, down, front, back = 6
    // +diags =8
    // 3D diags =27
    // calc center of mass tile
    // max distance is then L/2

    // only receiving needs to be unambiguous; in terms of tiles in node

    return 0;
  }


}; // end of Node class

} // end of corgi namespace


