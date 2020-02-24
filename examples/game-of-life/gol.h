#pragma once

#include <vector>

#include "../../tile.h"
#include "../../corgi.h"
#include "../../toolbox/dataContainer.h"

#include <mpi4cpp/mpi.h>


namespace gol {



/// Snapshot patch of a CA simulation
class Mesh {

  public:

  /// patch dimensions
  int Nx;
  int Ny;

  int halo = 1;

  /// Internal 2D mesh storing the values
  std::vector<int> mesh;

  /// ctor
  Mesh(int Nx, int Ny);

  // Indexing with +1 halo regions around the array
  int indx(const int i, const int j) const {
    return (i+halo) + (Nx +2*halo)*(j+halo);
  }

  int size()
  {
    return ((Nx + 2*halo) * (Ny + 2*halo));    
  }

  void clear()
  {
    std::fill(mesh.begin(), mesh.end(), 0);
  }

  /// 2D access operator for values
  int operator () (const int i, const int j) const {
    return mesh[ indx(i,j) ];
  }

  int &operator () (const int i, const int j) {
    return mesh[ indx(i,j) ];
  }


  void copy_vert(Mesh& rhs, int lhsI, int rhsI);

  void copy_horz(Mesh& rhs, int lhsJ, int rhsJ);


};


/// Simulation grid
//class Grid : public corgi::Grid {
//
//  public:
//    Grid(size_t nx, size_t ny) : corgi::Grid(nx, ny) { }
//
//    ~Grid() { };
//
//    // std::string pet_shop();
//
//    /// Cycle data containers of each cell forward
//    // void cycle() {
//    //   for (auto& it: cells) {
//    //     auto cellptr = std::dynamic_pointer_cast<CA_tile>( it.second );
//    //     cellptr->data.cycle();
//    //   }
//    // }
//
//};



/// Small local cellular automata patch
class Tile : public corgi::Tile<2> {

  public:

    using Tile_t = Tile;
    using Tileptr = std::shared_ptr<Tile>;

    Tile() = default;

    ~Tile() override = default;

    // extending the base class
    datarotators::Rotator<Mesh,2> data;

    void add_data(Mesh m);

    Mesh& get_data();

    Mesh* get_dataptr();

    Mesh& get_new_data();

    void update_boundaries(corgi::Grid<2>& grid);

    /// step forward
    void cycle() { data.cycle(); }

    std::vector<mpi4cpp::mpi::request> 
    send_data( mpi4cpp::mpi::communicator&, int orig, int mode, int tag) override;

    std::vector<mpi4cpp::mpi::request> 
    recv_data( mpi4cpp::mpi::communicator&, int dest, int mode, int tag) override;

};



class Solver {

  public:
    void solve(Tile&);

};





} // end of namespace gol
