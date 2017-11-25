#include <string>

#include "gol.h"

using namespace gol;

/// initialize internal mesh
patch::patch(size_t Nx, size_t Ny) : Nx(Nx), Ny(Ny) {
  mesh.resize(Nx * Ny);    
}


std::string CellularAutomataCell::bark() { return "Woof!"; };


// std::string Grid::petShop() { return "No Corgis for sale."; };





