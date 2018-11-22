#include <string>

#include "gol.h"
#include "../../toolbox/dataContainer.h"


using namespace gol;

/// initialize internal mesh
Mesh::Mesh(int Nx, int Ny) : Nx(Nx), Ny(Ny) {
  mesh.resize((Nx + 2*halo) * (Ny + 2*halo));    
}


/// Copy vertical slice
void Mesh::copy_vert(Mesh& rhs, int lhsI, int rhsI) {
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int j=0; j<this->Ny; j++) { 
    this->operator()(lhsI, j) = rhs(rhsI, j);
  }
}


/// Copy horizontal slice 
void Mesh::copy_horz(Mesh& rhs, int lhsJ, int rhsJ) {
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");

  for(int i=0; i<this->Nx; i++) { 
    this->operator()(i, lhsJ) = rhs(i, rhsJ);
  }
}


// --------------------------------------------------


/// Add data to the container
void CA_tile::add_data(Mesh m) {
  data.push_back(m);
}


/// get current patch
Mesh& CA_tile::get_data() {
  return *data.get();
};

Mesh* CA_tile::get_dataptr() {
  return data.get();
};

/// get new data
Mesh& CA_tile::get_new_data() {
  return *data.get_new();
};


/// Update boundary/halo regions from neighbors
void CA_tile::update_boundaries(Grid& grid) {

  Mesh& mesh = get_data(); // target as a reference to update into
  

  // left 
  Tileptr cleft = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(-1, 0) ));
  Mesh& mleft = cleft->get_data();
  mesh.copy_vert(mleft, -1, mleft.Nx-1); // copy from right side to left

  // right
  Tileptr cright = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(+1, 0) ));
  Mesh& mright = cright->get_data();
  mesh.copy_vert(mright, mesh.Nx, 0); // copy from left side to right

  // top 
  Tileptr ctop = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(0, +1) ));
  Mesh& mtop = ctop->get_data();
  mesh.copy_horz(mtop, mesh.Ny, 0); // copy from bottom side to top

  // bottom
  Tileptr cbot = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(0, -1) ));
  Mesh& mbot = cbot->get_data();
  mesh.copy_horz(mbot, -1, mbot.Ny-1); // copy from top side to bottom




  // diagonals/corners
  // --------------------------------------------------  
    
  // top right
  Tileptr ctr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( this->neighs(+1, +1) ));
  Mesh& mtr = ctr->get_data();
  mesh(mesh.Nx, mesh.Ny) = mtr(0,0);
  

  // top left
  Tileptr ctl = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( this->neighs(-1, +1) ));
  Mesh& mtl = ctl->get_data();
  mesh(-1, mesh.Ny) = mtl(mtl.Nx-1,0);


  // bottom right
  Tileptr cbr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( this->neighs(+1, -1) ));
  Mesh& mbr = cbr->get_data();
  mesh(mesh.Nx, -1) = mbr(0,mbr.Ny-1);

  // bottom left WORKS
  Tileptr cbl = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( this->neighs(-1, -1) ));
  Mesh& mbl = cbl->get_data();
  mesh(-1,-1) = mbl(mbl.Nx-1, mbl.Ny-1);
  

};


void Solver::solve(CA_tile& cell) {
  Mesh& m    = cell.get_data();
  Mesh& mnew = cell.get_new_data();


  for(int i=0; i<(int)m.Nx; i++) { 
    for(int j=0; j<(int)m.Ny; j++) { 

      // count alives
      int alive = 0;
      for(int ir=-1; ir<2; ir++) {
        for(int jr=-1; jr<2; jr++) {
          int state = m( i + ir, j + jr );

          if(state == 1) { alive++; }

        }
      }

      // apply rules
      if(alive == 3) {
        mnew(i,j) = 1;
      // } else if(alive != 2) {
      //   mnew(i,j) = 0;
      }


    }
  }


}



// std::string Grid::pet_shop() { return "No Corgis for sale."; };

