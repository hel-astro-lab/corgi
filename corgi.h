#pragma once

#include <fmt/format.h>
#include <fmt/format.cc>
#include <fmt/string.h>
#include <fmt/ostream.h>

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef> // for offsetof
#include <unordered_map>

#include "mpi.h"

// #include "definitions.h"
#include "common.h"
#include "SparseGrid.h"



namespace corgi {

    class Cell {

        public:
            // Order here is fixed for mpi_cell_t

            /// unique cell ID
            uint64_t cid;

            /// MPI rank of who owns me
            int owner;

            /// coarse mpiGrid grid indices
            size_t i, j;

            /// If I am a virtual cell, who do I share the values the most.
            int top_virtual_owner;

            /// how many times do I have to be sent to others
            size_t communications;

            /// How many virtual neighbors do I have
            size_t number_of_virtual_neighbors = 0;

            /// Cell type listing
            bool local;

            std::vector<int> types;


            /// initalize cell according to its location (i,j) and owner (o)
            Cell(size_t i, size_t j, int o) {
                this->i     = i;
                this->j     = j;
                this->owner = o;
            }

            /*! \brief *virtual* base class destructor 
             * NOTE: this needs to be virtual so that child classes can be 
             * destroyed.
             */
            virtual ~Cell() { }

            /// return mpiGrid index
            const std::tuple<size_t, size_t> index() {
                return std::make_tuple( i, j );
            }

            /// return index of cells in relative to my position
            const std::tuple<size_t, size_t> neighs(int ir, int jr) {
                size_t ii = BC::xwrap( (int)this->i + ir );
                size_t jj = BC::ywrap( (int)this->j + jr );
                return std::make_tuple( ii, jj );
            }


            /// Return full neighborhood around me
            std::vector< std::tuple<size_t, size_t> > nhood() {
                std::vector< std::tuple<size_t, size_t> > nh;
                for (int ir=-1; ir<=1; ir++) {
                    for (int jr=-1; jr<=1; jr++) {
                        if (!( ir == 0 && jr == 0 )) {
                            nh.push_back( neighs(ir, jr) );
                        }
                    }
                }
                return nh;
            }


            /// Check if cell fulfills a single criteria
            bool is_type( int criteria ) {
                if( std::find(
                            types.begin(), 
                            types.end(), 
                            criteria) 
                        == types.end() 
                  ) {
                    return false;
                } 
                return true;
            }
                
            /// Vectorized version requiring cell to fulfill every criteria
            bool is_types( std::vector<int> criteria ) {
                for (auto crit: criteria) {
                    if (is_type(crit))  {
                            continue;
                    } else {
                        return false;
                    }
                }

                // passed all criteria
                return true;
            }

    }; // end of Cell class



    class Node {


        /// Global large scale grid where information
        // of all the mpi processes is stored
        // int _mpiGrid[conf::Nx][conf::Ny];
        // std::vector<std::vector<int> > _mpiGrid;
        SparseGrid<int> _mpiGrid;


        public:

            /// Map with cellID & cell data
            std::unordered_map< uint64_t, std::shared_ptr<corgi::Cell> > cells;


            /// list of cell id's that are to be sent to others
            std::vector<uint64_t> send_queue;

            /// list containing lists to where the aforementioned send_queue cells are to be sent
            std::vector< std::vector<int> > send_queue_address;

            /// Return global grid sizes
            size_t getNx() { return conf::Nx; };
            size_t getNy() { return conf::Ny; };

            /// Return global grid dimensions
            double getXmin() { return conf::xmin; };
            double getXmax() { return conf::xmax; };
            double getYmin() { return conf::ymin; };
            double getYmax() { return conf::ymax; };


            /// get mpi process for whatever location
            int mpiGrid(const size_t i, const size_t j) {
                return _mpiGrid(i,j);
            }

            /// set new mpi process for some cell
            void setMpiGrid(const size_t i, const size_t j, int val) {
                _mpiGrid(i,j) = val;
            }


            /// Create unique cell ids based on Morton z ordering
            uint64_t cellId(size_t i, size_t j) {
                return uint64_t( j*conf::Nx + i );
            }
            
            /*
            uint64_t cell_id( std::tuple<size_t, size_t> indx ) {
                size_t i = std::get<0>(indx);
                size_t j = std::get<1>(indx);
                return uint64_t(i * conf::Nx) + uint64_t(j);
            }
            */

            /// Add local cell to the node
            void addLocalCell( corgi::Cell c ) {

                // calculate unique global cell ID
                uint64_t cid = cellId(c.i, c.j);

                //TODO Catch error if cell is not already mine?
                c.cid   = cid;
                c.owner = rank;
                c.local = true;
                // c.types.push_back( cellType::LOCAL );
                
                cells.insert( std::make_pair(cid, &c) );
            }


            /// Return pointer to the actual cell data
            corgi::Cell* getCellData(const uint64_t cid) const {
            	if (this->cells.count(cid) > 0) {
            		return (corgi::Cell*) &(this->cells.at(cid));
            	} else {
            		return NULL;
            	}
            }

            /// Same as get_cell_data but with additional syntax sugar 
            corgi::Cell* operator [] (const uint64_t cid) const {
                return getCellData(cid);
            }

            /// Get a *copy* of the full cell; this is not what one usually wants
            corgi::Cell getCell( uint64_t cid ) {
                return *cells.at(cid);
            }


            /*! Return a vector of cell indices that fulfill a given criteria.  */
            std::vector<uint64_t> getAllCells(
                    const std::vector<int>& criteria = std::vector<int>(),
                    const bool sorted=false ) {
                std::vector<uint64_t> ret;

                for (auto it: cells) {
                    if (criteria.size() == 0) {
                        ret.push_back( it.first );
                        continue;
                    }

                    // criteria checking
                    auto c = it.second;
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


            /// Return all local cells
            std::vector<uint64_t> getCells(
                    const std::vector<int>& criteria = std::vector<int>(),
                    const bool sorted=false ) {

                std::vector<uint64_t> cell_list = getAllCells(criteria, sorted);

                size_t i = 0, len = cell_list.size();
                while (i < len) {
                    if (!cells.at( cell_list[i] )->local) {
                        std::swap(cell_list[i], cell_list.back());
                        cell_list.pop_back();
                        len -= 1;
                    } else {
                        i++;
                    }
                }

                return cell_list;
            }


            /// Return all cells that are of VIRTUAL type.
            std::vector<uint64_t> getVirtuals(
                    const std::vector<int>& criteria = std::vector<int>(),
                    const bool sorted=false ) {
                std::vector<uint64_t> cell_list = getAllCells(criteria, sorted);

                size_t i = 0, len = cell_list.size();
                while (i < len) {
                    if (cells.at( cell_list[i] )->local) {
                        std::swap(cell_list[i], cell_list.back());
                        cell_list.pop_back();
                        len -= 1;
                    } else {
                        i++;
                    }
                }

                return cell_list;
            }
            

            /// Check if we have a cell with the given index
            // bool is_local(std::tuple<int, int> indx) {
            bool isLocal(uint64_t cid) {
                bool local = false;

                // Do we have it on the list=
                if (cells.count( cid ) > 0) {
                    // is it local (i.e., not virtual)
                    if ( cells.at(cid)->local ) {
                        local = true;
                    }
                }

                return local;
            }

            // TODO: relative indexing w.r.t. given cell
            // std::tuple<size_t, size_t> get_neighbor_index(corgi::Cell, int i, int j) {
            //     return c.neighs( std::make_tuple(i,j) );
            // }

            // TODO: get_neighbor_cell(c, i, j)



            std::vector<int> virtualNeighborhood(uint64_t cid) {

                auto c = getCellData(cid);
                std::vector< std::tuple<size_t, size_t> > neigs = c->nhood();
                std::vector<int> virtual_owners;
                for (auto indx: neigs) {

                    /* TODO: check boundary cells here; 
                     * now we assume periodicity in x and y
                    if (std::get<0>(indx) == ERROR_INDEX ||
                        std::get<1>(indx) == ERROR_INDEX) {
                        continue;
                    }
                    */

                    // Get cell id from index notation
                    size_t i = std::get<0>(indx);
                    size_t j = std::get<1>(indx);
                    uint64_t cid = cellId(i, j);

                    if (!isLocal( cid )) {
                        int whoami = _mpiGrid(indx); 
                        virtual_owners.push_back( whoami );
                    }
                }

                return virtual_owners;
            }
              
              
            // Number of virtual neighbors that the cell might have.
            /*
            size_t number_of_virtual_neighborhood(corgi::Cell c) {
                return virtual_neighborhood(c).size();
            }
            */


            /*! Analyze my local boundary cells that will be later on
             * send to the neighbors as virtual cells. 
             *
             * This is where the magic happens and we analyze what and who to send to.
             * These values *must* be same for everybody, this is why we use
             * mode of the owner list and in case of conflict pick the smaller value.
             * This way everybody knows what to expect and we avoid creating conflicts 
             * in communication. This information is then being sent to other processes 
             * together with the cells and is analyzed there by others inside the
             * `rank_virtuals` function.
             * */
            void analyzeBoundaryCells() {

                for (auto cid: getCells()) {
                    std::vector<int> virtual_owners = virtualNeighborhood(cid);
                    size_t N = virtual_owners.size();

                    // If N > 0 then this is a boundary cell.
                    // other criteria could also apply but here we assume
                    // neighborhood according to spatial distance.
                    if (N > 0) {

                        /* Now we analyze `owner` vector as:
                         * - sort the vector
                         * - compute mode of the list to see who owns most of the
                         * - remove repeating elements creating a unique list. */
                         
                        // sort
                        std::sort( virtual_owners.begin(), virtual_owners.end() );

                        // compute mode by creating a frequency array
                        // NOTE: in case of same frequency we implicitly pick smaller rank
                        int max=0, top_owner = virtual_owners[0];
                        for(size_t i=0; i<virtual_owners.size(); i++) {
                            int co = (int)count(virtual_owners.begin(), 
                                            virtual_owners.end(), 
                                            virtual_owners[i]);
                            if(co > max) {      
                                max = co;
                                top_owner = virtual_owners[i];
                            }
                        } 

                        // remove duplicates
                        virtual_owners.erase( unique( virtual_owners.begin(), 
                                              virtual_owners.end() 
                                            ), virtual_owners.end() );


                        // update cell values
                        auto c = getCellData(cid);
                        c->top_virtual_owner = top_owner;
                        c->communications    = virtual_owners.size();
                        c->number_of_virtual_neighbors = N;

                        if (std::find( send_queue.begin(),
                                       send_queue.end(),
                                       cid) == send_queue.end()
                           ) {
                            send_queue.push_back( cid );
                            send_queue_address.push_back( virtual_owners );
                        }
                    }
                }
            }
            

            /// Clear send queue, issue this only after the send has been successfully done
            void clearSendQueue() {
                send_queue.clear();
                send_queue_address.clear();
            }


        public:
            // -------------------------------------------------- 
            /// MPI communication related stuff
            int rank  = 0;
            int Nrank = 0;
            MPI_Comm comm;

            //TODO double definition for python debugging
            bool master = false;

            MPI_Datatype mpi_cell_t;

            std::vector<MPI_Request> sent_info_messages;
            std::vector<MPI_Request> sent_cell_messages;

            std::vector<MPI_Request> recv_info_messages;
            std::vector<MPI_Request> recv_cell_messages;

            Node() {
                fmt::print("initializing node ({} {})...\n", conf::Nx, conf::Ny);
                
                // allocating _mpiGrid
                _mpiGrid.resize(conf::Nx, conf::Ny);

            }


            /// Initialize MPI and related auxiliary variables
            void initMpi() {

                //--------------------------------------------------
                // Start MPI
                //
                // TODO do this in main program with arg and argv
                MPI_Init(NULL, NULL);

                comm = MPI_COMM_WORLD;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &Nrank);

                // detect master
                if (rank == MASTER_RANK) { master = true; };

                fmt::print("Hi from rank {}\n", rank);
                if (master) { fmt::print("master is {}\n", rank); };


                //--------------------------------------------------
                // Initialize the cell frame type
                int count = 7;
                int blocklens[] = { 1, 1, 1, 1, 1, 1, 1 };
                MPI_Aint indices[7];
                indices[0] = (MPI_Aint)offsetof( Cell, cid);
                indices[1] = (MPI_Aint)offsetof( Cell, owner);
                indices[2] = (MPI_Aint)offsetof( Cell, i);
                indices[3] = (MPI_Aint)offsetof( Cell, j);
                indices[4] = (MPI_Aint)offsetof( Cell, top_virtual_owner);
                indices[5] = (MPI_Aint)offsetof( Cell, communications);
                indices[6] = (MPI_Aint)offsetof( Cell, number_of_virtual_neighbors);
                
                MPI_Datatype types[] = {
                                      MPI_UINT64_T,  // cid
                                      MPI_INT,       // owner
                                      MPI_SIZE_T,    // i
                                      MPI_SIZE_T,    // j
                                      MPI_INT,       // top_virtual_owner
                                      MPI_SIZE_T,    // communications
                                      MPI_SIZE_T     // num. of virt. owners.
                                       };
                MPI_Type_create_struct(count, blocklens, indices, types, &mpi_cell_t);
                MPI_Type_commit(&mpi_cell_t);
                
                //--------------------------------------------------

            }


            /// Finalize MPI environment 
            void finalizeMpi() {
                MPI_Type_free(&mpi_cell_t);
                MPI_Finalize();
            }


            /// Broadcast master ranks mpiGrid to everybody
            void bcastMpiGrid() {

                std::vector<int> tmp;
                if (master) {
                  tmp = _mpiGrid.serialize();
                } else {
                  tmp.resize(conf::Nx * conf::Ny);
                  for(size_t k=0; k<conf::Nx*conf::Ny; k++) {tmp[k] = -1.0;};
                }


                MPI_Bcast(&tmp[0],
                          conf::Nx*conf::Ny, 
                          MPI_INT, 
                          MASTER_RANK, 
                          MPI_COMM_WORLD
                         );

                // unpack
                if(!master) {
                  _mpiGrid.unpack(tmp, conf::Nx, conf::Ny);
                }

            }

            /// Issue isends to everywhere
            // First we send a warning message of how many cells to expect.
            // Based on this the receiving side can prepare accordingly.
            void communicateSendCells() {

                sent_info_messages.clear();
                sent_cell_messages.clear();
                int j = 0;
                
                for (int dest = 0; dest<Nrank; dest++) {
                    if(dest == rank) { continue; } // do not send to myself
                    
                    int i = 0;
                    std::vector<int> to_be_sent;
                    for (std::vector<int> address: send_queue_address) {
                        if( std::find( address.begin(),
                                       address.end(),
                                       dest) != address.end()) {
                            to_be_sent.push_back( i );
                        }
                        i++;
                    }

                    // initial message informing how many cells are coming
                    // TODO: this whole thing could be avoided by using 
                    // MPI_Iprobe in the receiving end. Maybe.
                    uint64_t Nincoming_cells = uint64_t(to_be_sent.size());

                    MPI_Request req;
                    sent_info_messages.push_back( req );

                    MPI_Isend(
                            &Nincoming_cells, 
                            1,
                            MPI_UNSIGNED_LONG_LONG,
                            dest,
                            commType::NCELLS,
                            comm,
                            &sent_info_messages[j] 
                            );
                    j++;
                }


                // send the real cell data now
                // We optimize this by only packing the cell data
                // once, and then sending the same thing to everybody who needs it.
                int i = 0;
                for (auto cid: send_queue) {
                    sendCellData( cid, send_queue_address[i] );
                    i++;
                }

            }


            /// Pack cell and send to everybody on the dests list
            void sendCellData(uint64_t cid, std::vector<int> dests) {
                auto c = getCellData(cid);
                
                size_t j = sent_cell_messages.size();
                
                for (auto dest: dests) {
                    MPI_Request req;
                    sent_cell_messages.push_back( req );

                    MPI_Isend(
                            c,
                            1,
                            mpi_cell_t,
                            dest,
                            commType::CELLDATA,
                            comm,
                            &sent_cell_messages[j]
                            );
                    j++;
                }
            }


            /// Receive incoming stuff
            void communicateRecvCells() {

                recv_info_messages.clear();
                recv_cell_messages.clear();

                size_t i = 0;
                for (int source=0; source<Nrank; source++) {
                    if (source == rank) { continue; } // do not receive from myself

                    // communicate with how many cells there are incoming

                    // TODO: use MPI_IProbe to check if there are 
                    // any messages for me instead of assuming that there is

                    MPI_Request req;
                    recv_info_messages.push_back( req );

                    uint64_t Nincoming_cells;
                    MPI_Irecv(
                            &Nincoming_cells,
                            1,
                            MPI_UNSIGNED_LONG_LONG,
                            source,
                            commType::NCELLS,
                            comm,
                            &recv_info_messages[i]
                            );

                    // TODO: Remove this code block and do in background instead
                    MPI_Wait(&recv_info_messages[i], MPI_STATUS_IGNORE);
                    
                    /*
                    fmt::print("{}: I got a message! Waiting {} cells from {}\n",
                            rank, Nincoming_cells, source);
                    */


                    // Now receive the cells themselves
                    size_t j = recv_cell_messages.size();
                    for (size_t ic=0; ic<Nincoming_cells; ic++) {
                        Cell inc_c(0,0,0); // TODO: initialize with better default values

                        MPI_Request reqc;
                        recv_cell_messages.push_back( reqc );
                        MPI_Irecv(
                                &inc_c,
                                1,
                                mpi_cell_t,
                                source,
                                commType::CELLDATA,
                                comm,
                                &recv_cell_messages[j]
                                );

                        MPI_Wait(&recv_cell_messages[j], MPI_STATUS_IGNORE);
                        j++;

                        uint64_t cid = inc_c.cid;
            	        if (this->cells.count(cid) == 0) {
                            // Cell does not exist yet; create it
                            // TODO: Check validity of the cell better
                            // TODO: Use add_cell() instead of directly 
                            //       probing the class interiors
                            
                            inc_c.local = false;
                            cells.insert( std::make_pair(cid, &inc_c) );

                        } else {
                            // Cell is already on my virtual list; update
                            // TODO: use = operator instead.
                            auto c = getCellData(cid);

                            if (c->local) {
                                // TODO: better error handling; i.e. resolve the conflict
                                fmt::print("{}: ERROR trying to add virtual cell that is already local\n", rank);
                                exit(1);
                            }

                            c->owner             = inc_c.owner;
                            c->i                 = inc_c.i;
                            c->j                 = inc_c.j;
                            c->top_virtual_owner = inc_c.top_virtual_owner;
                            c->communications    = inc_c.communications;
                            c->number_of_virtual_neighbors = inc_c.number_of_virtual_neighbors;


                        };

                    }
                    i++;
                }
            }

    }; // end of Node class

} // end of corgi namespace


