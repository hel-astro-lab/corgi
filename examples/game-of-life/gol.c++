#include <string>

#include "gol.h"

using namespace gol;

/// initialize internal mesh
Patch::Patch(size_t Nx, size_t Ny) : Nx(Nx), Ny(Ny) {
  mesh.resize(Nx * Ny);    
}




// std::string Grid::petShop() { return "No Corgis for sale."; };





