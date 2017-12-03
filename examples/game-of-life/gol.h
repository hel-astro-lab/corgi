#pragma once

#include <vector>

#include "../../cell.h"
#include "../../corgi.h"
#include "../../toolbox/dataContainer.h"



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

  /// 2D access operator for values
  int operator () (const int i, const int j) const {
    return mesh[ indx(i,j) ];
  }

  int &operator () (const int i, const int j) {
    return mesh[ indx(i,j) ];
  }


  void copyVert(Mesh& rhs, int lhsI, int rhsI);

  void copyHorz(Mesh& rhs, int lhsJ, int rhsJ);


};


/// Simulation grid
class Grid : public corgi::Node {

  public:
    Grid(size_t nx, size_t ny) : corgi::Node(nx, ny) { }

    ~Grid() { };

    // std::string petShop();

    /// Cycle data containers of each cell forward
    // void cycle() {
    //   for (auto& it: cells) {
    //     auto cellptr = std::dynamic_pointer_cast<CellularAutomataCell>( it.second );
    //     cellptr->data.cycle();
    //   }
    // }

};



/// Small local cellular automata patch
class CellularAutomataCell : public corgi::Cell {

  public:

    typedef CellularAutomataCell CellType;
    typedef std::shared_ptr<CellularAutomataCell> CellPtr;

    CellularAutomataCell(size_t i, size_t j, 
             int o, 
             size_t nx, size_t ny
             ) : corgi::Cell(i, j, o, nx, ny) { }

    ~CellularAutomataCell() { };

    // extending the base class
    datarotators::DataContainer<Mesh> data;

    void addData(Mesh m);

    Mesh& getData();

    Mesh* getDataPtr();

    Mesh& getNewData();

    void updateBoundaries(Grid& grid);


    /// step forward
    void cycle() { data.cycle(); }

};



class Solver {


  public:
    void solve(CellularAutomataCell& cell);



};





} // end of namespace gol
