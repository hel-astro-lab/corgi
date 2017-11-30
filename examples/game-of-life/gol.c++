#include <string>

#include "gol.h"

using namespace gol;

/// initialize internal mesh
Mesh::Mesh(size_t Nx, size_t Ny) : Nx(Nx), Ny(Ny) {
  mesh.resize(Nx * Ny);    
}





// std::string Grid::petShop() { return "No Corgis for sale."; };

