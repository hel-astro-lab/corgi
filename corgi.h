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

// #include <cstddef> // for offsetof

#include "mpi.h"

// #include "common.h"
#include "toolbox/SparseGrid.h"
#include "cell.h"


namespace corgi {

class Node {

  private:

  // --------------------------------------------------
  /// Global grid dimensions
  size_t Nx = 0;
  size_t Ny = 0;

  /// Global simulation box size
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;

  public:

  /// Return global grid sizes
  size_t getNx() { return Nx; };
  size_t getNy() { return Ny; };

  /// Return global grid dimensions
  double getXmin() { return xmin; };
  double getXmax() { return xmax; };
  double getYmin() { return ymin; };
  double getYmax() { return ymax; };


  /// Set physical grid size
  void setGridLims(double _xmin, double _xmax, 
                   double _ymin, double _ymax) {
    xmin = _xmin;
    xmax = _xmax;

    ymin = _ymin;
    ymax = _ymax;
  }


  private:
  // --------------------------------------------------
  // Cell Mapping
  typedef corgi::Cell                           CellType;
  // typedef std::unique_ptr<CellType>             CellPtr;
  typedef std::shared_ptr<CellType>             CellPtr;
  typedef std::unordered_map<uint64_t, CellPtr> CellMap;

  public:
  /// Map with cellID & cell data
  CellMap cells;

  /*! Global large scale block grid where information
   * of all the mpi processes are stored
   */
  SparseGrid<int> mpiGrid;

  // Python bindings for mpiGrid
  int pyGetMpiGrid(size_t i, size_t j) {
    return mpiGrid(i,j);
  }

  void pySetMpiGrid(size_t i, size_t j, int val) {
    mpiGrid(i,j) = val;
  }

  /// Create unique cell ids based on Morton z ordering
  uint64_t cellId(size_t i, size_t j) {
    return uint64_t( j*Nx + i );
  }


  public:
  // --------------------------------------------------
  // Cell constructors & destructors
    
  /// Constructor
  Node(size_t nx, size_t ny) : Nx(nx), Ny(ny) {
    // fmt::print("initializing node ({} {})...\n", Nx, Ny);
      
    // allocating _mpiGrid
    mpiGrid.resize(Nx, Ny);
  }

  /// Deallocate and free everything
  ~Node() { } 


  public:
  // --------------------------------------------------
  // Cell addition etc. manipulation
    
  /// Add local cell to the node
  // void addCell(CellType& cell) {
  void addCell(CellPtr cellptr) {

    // claim unique ownership of the cell (for unique_ptr)
    // std::unique_ptr<corgi::Cell> cellptr = std::make_unique<corgi:Cell>(cell);
    // CellPtr cellptr = std::make_unique<CellType>(cell);
    
    // calculate unique global cell ID
    uint64_t cid = cellId(cellptr->i, cellptr->j);
  
    cellptr->cid   = cid;
    cellptr->owner = rank;
    cellptr->local = true; //TODO Catch error if cell is not already mine?

    // cells.emplace(cid, std::move(cellptr)); // unique_ptr needs to be moved
    cells.emplace(cid, cellptr); // NOTE using c++14 emplace to avoid copying
  }


  /*! Return a vector of cell indices that fulfill a given criteria.  */
  std::vector<uint64_t> getCellIds(
      const std::vector<int>& criteria = std::vector<int>(),
      const bool sorted=false ) {
    std::vector<uint64_t> ret;

    for (auto& it: cells) {
      if (criteria.size() == 0) {
        ret.push_back( it.first );
        continue;
      }

      // criteria checking
      auto& c = it.second;
      if (!c->is_types( criteria ) ) {
        continue;
      }

      ret.push_back( it.first );
    }


    // optional sort based on the cell id
    if (sorted && ret.size() > 0) {
      std::sort(ret.begin(), ret.end());
    }

    return ret;
  }

  /*! \brief Get individual cell (as a reference)
   *
   * NOTE: from StackOverflow (recommended getter method):
   * OtherType& get_othertype(const std::string& name)
   * {
   *     auto it = otMap.find(name);
   *     if (it == otMap.end()) throw std::invalid_argument("entry not found");
   *     return *(it->second);
   * }
   *
   * This way map retains its ownership of the cell and we avoid giving pointers
   * away from the Class.
   */
  CellType& getCell(const uint64_t cid) {
    auto it = cells.find(cid);
    if (it == cells.end()) throw std::invalid_argument("cell entry not found");

    return *(it->second);
  }

  /// \brief Get individual cell (as a pointer)
  CellPtr getCellPtr(const uint64_t cid) {
    auto it = cells.find(cid);
    if (it == cells.end()) throw std::invalid_argument("entry not found");
    return it->second;
  }




  // /// Return pointer to the actual cell data
  // corgi::Cell* getCellData(const uint64_t cid) const {
  //   if (this->cells.count(cid) > 0) {
  //     return (corgi::Cell*) &(this->cells.at(cid));
  //   } else {
  //     return NULL;
  //   }
  // }

  // /// Same as get_cell_data but with additional syntax sugar 
  // corgi::Cell* operator [] (const uint64_t cid) const {
  //   return getCellData(cid);
  // }

  // /// Get a *copy* of the full cell; this is not what one usually wants
  // corgi::Cell getCell( uint64_t cid ) {
  //   return *cells.at(cid);
  // }

  // /// Return all local cells
  // std::vector<uint64_t> getCells(
  //     const std::vector<int>& criteria = std::vector<int>(),
  //     const bool sorted=false ) {

  //   std::vector<uint64_t> cell_list = getAllCells(criteria, sorted);

  //   size_t i = 0, len = cell_list.size();
  //   while (i < len) {
  //     if (!cells.at( cell_list[i] )->local) {
  //       std::swap(cell_list[i], cell_list.back());
  //       cell_list.pop_back();
  //       len -= 1;
  //     } else {
  //       i++;
  //     }
  //   }

  //   return cell_list;
  // }


  // /// Return all cells that are of VIRTUAL type.
  // std::vector<uint64_t> getVirtuals(
  //     const std::vector<int>& criteria = std::vector<int>(),
  //     const bool sorted=false ) {
  //   std::vector<uint64_t> cell_list = getAllCells(criteria, sorted);

  //   size_t i = 0, len = cell_list.size();
  //   while (i < len) {
  //     if (cells.at( cell_list[i] )->local) {
  //       std::swap(cell_list[i], cell_list.back());
  //       cell_list.pop_back();
  //       len -= 1;
  //     } else {
  //       i++;
  //     }
  //   }

  //   return cell_list;
  // }


  // /// Check if we have a cell with the given index
  // // bool is_local(std::tuple<int, int> indx) {
  // bool isLocal(uint64_t cid) {
  //   bool local = false;

  //   // Do we have it on the list=
  //   if (cells.count( cid ) > 0) {
  //     // is it local (i.e., not virtual)
  //     if ( cells.at(cid)->local ) {
  //       local = true;
  //     }
  //   }

  //   return local;
  // }

  // // TODO: relative indexing w.r.t. given cell
  // // std::tuple<size_t, size_t> get_neighbor_index(corgi::Cell, int i, int j) {
  // //     return c.neighs( std::make_tuple(i,j) );
  // // }

  // // TODO: get_neighbor_cell(c, i, j)



  // std::vector<int> virtualNeighborhood(uint64_t cid) {

  //   auto c = getCellData(cid);
  //   std::vector< std::tuple<size_t, size_t> > neigs = c->nhood();
  //   std::vector<int> virtual_owners;
  //   for (auto indx: neigs) {

  //     /* TODO: check boundary cells here; 
  //      * now we assume periodicity in x and y
  //      if (std::get<0>(indx) == ERROR_INDEX ||
  //      std::get<1>(indx) == ERROR_INDEX) {
  //      continue;
  //      }
  //      */

  //     // Get cell id from index notation
  //     size_t i = std::get<0>(indx);
  //     size_t j = std::get<1>(indx);
  //     uint64_t cid = cellId(i, j);

  //     if (!isLocal( cid )) {
  //       int whoami = _mpiGrid(indx); 
  //       virtual_owners.push_back( whoami );
  //     }
  //   }

  //   return virtual_owners;
  // }


  // // Number of virtual neighbors that the cell might have.
  // /*
  //    size_t number_of_virtual_neighborhood(corgi::Cell c) {
  //    return virtual_neighborhood(c).size();
  //    }
  //    */


  // /*! Analyze my local boundary cells that will be later on
  //  * send to the neighbors as virtual cells. 
  //  *
  //  * This is where the magic happens and we analyze what and who to send to.
  //  * These values *must* be same for everybody, this is why we use
  //  * mode of the owner list and in case of conflict pick the smaller value.
  //  * This way everybody knows what to expect and we avoid creating conflicts 
  //  * in communication. This information is then being sent to other processes 
  //  * together with the cells and is analyzed there by others inside the
  //  * `rank_virtuals` function.
  //  * */
  // void analyzeBoundaryCells() {

  //   for (auto cid: getCells()) {
  //     std::vector<int> virtual_owners = virtualNeighborhood(cid);
  //     size_t N = virtual_owners.size();

  //     // If N > 0 then this is a boundary cell.
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


  //       // update cell values
  //       auto c = getCellData(cid);
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
  // void clearSendQueue() {
  //   send_queue.clear();
  //   send_queue_address.clear();
  // }



  // // --------------------------------------------------
  // // Send queues etc.
  //   
  // /// list of cell id's that are to be sent to others
  // std::vector<uint64_t> send_queue;

  // /// list containing lists to where the aforementioned send_queue cells are to be sent
  // std::vector< std::vector<int> > send_queue_address;



  // public:
  // // -------------------------------------------------- 
  // /// MPI communication related stuff
  int rank  = 0;
  int Nrank = 0;
  MPI_Comm comm;

  // indicate master node
  bool master = false;

  // MPI_Datatype mpi_cell_t;

  // std::vector<MPI_Request> sent_info_messages;
  // std::vector<MPI_Request> sent_cell_messages;

  // std::vector<MPI_Request> recv_info_messages;
  // std::vector<MPI_Request> recv_cell_messages;





  // /// Initialize MPI and related auxiliary variables
  // void initMpi() {

  //   //--------------------------------------------------
  //   // Start MPI
  //   //
  //   // TODO do this in main program with arg and argv
  //   MPI_Init(NULL, NULL);

  //   comm = MPI_COMM_WORLD;
  //   MPI_Comm_rank(comm, &rank);
  //   MPI_Comm_size(comm, &Nrank);

  //   // detect master
  //   if (rank == MASTER_RANK) { master = true; };

  //   // fmt::print("Hi from rank {}\n", rank);
  //   // if (master) { fmt::print("master is {}\n", rank); };


  //   //--------------------------------------------------
  //   // Initialize the cell frame type
  //   int count = 7;
  //   int blocklens[] = { 1, 1, 1, 1, 1, 1, 1 };
  //   MPI_Aint indices[7];
  //   indices[0] = (MPI_Aint)offsetof( Cell, cid);
  //   indices[1] = (MPI_Aint)offsetof( Cell, owner);
  //   indices[2] = (MPI_Aint)offsetof( Cell, i);
  //   indices[3] = (MPI_Aint)offsetof( Cell, j);
  //   indices[4] = (MPI_Aint)offsetof( Cell, top_virtual_owner);
  //   indices[5] = (MPI_Aint)offsetof( Cell, communications);
  //   indices[6] = (MPI_Aint)offsetof( Cell, number_of_virtual_neighbors);

  //   MPI_Datatype types[] = {
  //     MPI_UINT64_T,  // cid
  //     MPI_INT,       // owner
  //     MPI_SIZE_T,    // i
  //     MPI_SIZE_T,    // j
  //     MPI_INT,       // top_virtual_owner
  //     MPI_SIZE_T,    // communications
  //     MPI_SIZE_T     // num. of virt. owners.
  //   };
  //   MPI_Type_create_struct(count, blocklens, indices, types, &mpi_cell_t);
  //   MPI_Type_commit(&mpi_cell_t);

  //   //--------------------------------------------------

  // }


  // /// Finalize MPI environment 
  // void finalizeMpi() {
  //   MPI_Type_free(&mpi_cell_t);
  //   MPI_Finalize();
  // }


  // /// Broadcast master ranks mpiGrid to everybody
  // void bcastMpiGrid() {

  //   std::vector<int> tmp;
  //   if (master) {
  //     tmp = _mpiGrid.serialize();
  //   } else {
  //     tmp.resize(Nx * Ny);
  //     for(size_t k=0; k<Nx*Ny; k++) {tmp[k] = -1.0;};
  //   }


  //   MPI_Bcast(&tmp[0],
  //       Nx*Ny, 
  //       MPI_INT, 
  //       MASTER_RANK, 
  //       MPI_COMM_WORLD
  //       );

  //   // unpack
  //   if(!master) {
  //     _mpiGrid.unpack(tmp, Nx, Ny);
  //   }

  // }

  // /// Issue isends to everywhere
  // // First we send a warning message of how many cells to expect.
  // // Based on this the receiving side can prepare accordingly.
  // void communicateSendCells() {

  //   sent_info_messages.clear();
  //   sent_cell_messages.clear();
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

  //     // initial message informing how many cells are coming
  //     // TODO: this whole thing could be avoided by using 
  //     // MPI_Iprobe in the receiving end. Maybe.
  //     uint64_t Nincoming_cells = uint64_t(to_be_sent.size());

  //     MPI_Request req;
  //     sent_info_messages.push_back( req );

  //     MPI_Isend(
  //         &Nincoming_cells, 
  //         1,
  //         MPI_UNSIGNED_LONG_LONG,
  //         dest,
  //         commType::NCELLS,
  //         comm,
  //         &sent_info_messages[j] 
  //         );
  //     j++;
  //   }


  //   // send the real cell data now
  //   // We optimize this by only packing the cell data
  //   // once, and then sending the same thing to everybody who needs it.
  //   int i = 0;
  //   for (auto cid: send_queue) {
  //     sendCellData( cid, send_queue_address[i] );
  //     i++;
  //   }

  // }


  // /// Pack cell and send to everybody on the dests list
  // void sendCellData(uint64_t cid, std::vector<int> dests) {
  //   auto c = getCellData(cid);

  //   size_t j = sent_cell_messages.size();

  //   for (auto dest: dests) {
  //     MPI_Request req;
  //     sent_cell_messages.push_back( req );

  //     MPI_Isend(
  //         c,
  //         1,
  //         mpi_cell_t,
  //         dest,
  //         commType::CELLDATA,
  //         comm,
  //         &sent_cell_messages[j]
  //         );
  //     j++;
  //   }
  // }


  // /// Receive incoming stuff
  // void communicateRecvCells() {

  //   recv_info_messages.clear();
  //   recv_cell_messages.clear();

  //   size_t i = 0;
  //   for (int source=0; source<Nrank; source++) {
  //     if (source == rank) { continue; } // do not receive from myself

  //     // communicate with how many cells there are incoming

  //     // TODO: use MPI_IProbe to check if there are 
  //     // any messages for me instead of assuming that there is

  //     MPI_Request req;
  //     recv_info_messages.push_back( req );

  //     uint64_t Nincoming_cells;
  //     MPI_Irecv(
  //         &Nincoming_cells,
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
  //        fmt::print("{}: I got a message! Waiting {} cells from {}\n",
  //        rank, Nincoming_cells, source);
  //        */


  //     // Now receive the cells themselves
  //     size_t j = recv_cell_messages.size();
  //     for (size_t ic=0; ic<Nincoming_cells; ic++) {
  //       Cell inc_c(0,0,0, Nx, Ny); // TODO: initialize with better default values

  //       MPI_Request reqc;
  //       recv_cell_messages.push_back( reqc );
  //       MPI_Irecv(
  //           &inc_c,
  //           1,
  //           mpi_cell_t,
  //           source,
  //           commType::CELLDATA,
  //           comm,
  //           &recv_cell_messages[j]
  //           );

  //       MPI_Wait(&recv_cell_messages[j], MPI_STATUS_IGNORE);
  //       j++;

  //       uint64_t cid = inc_c.cid;
  //       if (this->cells.count(cid) == 0) {
  //         // Cell does not exist yet; create it
  //         // TODO: Check validity of the cell better
  //         // TODO: Use add_cell() instead of directly 
  //         //       probing the class interiors

  //         inc_c.local = false;
  //         cells.insert( std::make_pair(cid, &inc_c) );

  //       } else {
  //         // Cell is already on my virtual list; update
  //         // TODO: use = operator instead.
  //         auto c = getCellData(cid);

  //         if (c->local) {
  //           // TODO: better error handling; i.e. resolve the conflict
  //           // TODO: use throw exceptions
  //           // fmt::print("{}: ERROR trying to add virtual cell that is already local\n", rank);
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


