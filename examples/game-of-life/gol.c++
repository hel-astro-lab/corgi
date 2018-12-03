#include <string>

#include "gol.h"
#include "../../toolbox/dataContainer.h"


using namespace gol;
using namespace mpi4cpp;

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
void Tile::add_data(Mesh m) {
  data.push_back(m);
}


/// get current patch
Mesh& Tile::get_data() {
  return data.get();
};

/// get new data
Mesh& Tile::get_new_data() {
  return data.get(1);
};


/// Update boundary/halo regions from neighbors
void Tile::update_boundaries(corgi::Node<2>& grid) {

  Mesh& mesh = get_data(); // target as a reference to update into

  // left 
  Tileptr cleft = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr(  neighs(-1, 0) )); 
  if(cleft) {
    Mesh& mleft = cleft->get_data();
    mesh.copy_vert(mleft, -1, mleft.Nx-1); // copy from right side to left
  }

  // right
  Tileptr cright = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(+1, 0) ));
  if(cright) {
    Mesh& mright = cright->get_data();
    mesh.copy_vert(mright, mesh.Nx, 0); // copy from left side to right
  }

  // top 
  Tileptr ctop = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr(   neighs(0, +1) ));
  if(ctop) {
    Mesh& mtop = ctop->get_data();
    mesh.copy_horz(mtop, mesh.Ny, 0); // copy from bottom side to top
  }

  // bottom
  Tileptr cbot = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr(   neighs(0, -1) ));
  if(cbot) {
    Mesh& mbot = cbot->get_data();
    mesh.copy_horz(mbot, -1, mbot.Ny-1); // copy from top side to bottom
  }



  // diagonals/corners
  // --------------------------------------------------  
    
  // top right
  Tileptr ctr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(+1, +1) ));
  if(ctr) {
    Mesh& mtr = ctr->get_data();
    mesh(mesh.Nx, mesh.Ny) = mtr(0,0);
  }

  // top left
  Tileptr ctl = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(-1, +1) ));
  if(ctl) {
    Mesh& mtl = ctl->get_data();
    mesh(-1, mesh.Ny) = mtl(mtl.Nx-1,0);
  }

  // bottom right
  Tileptr cbr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(+1, -1) ));
  if(cbr) {
    Mesh& mbr = cbr->get_data();
    mesh(mesh.Nx, -1) = mbr(0,mbr.Ny-1);
  }

  // bottom left
  Tileptr cbl = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(-1, -1) ));
  if(cbl) {
    Mesh& mbl = cbl->get_data();
    mesh(-1,-1) = mbl(mbl.Nx-1, mbl.Ny-1);
  }
  

};


mpi::request Tile::send_data( mpi::communicator& comm, int dest, int tag)
{
  //std::cout << "SEND to " << dest << "\n";
  Mesh& mesh = get_data(); 

  mpi::request req;
  req = comm.isend(dest, cid, mesh.mesh.data(), mesh.size() );

  return req;
}

mpi::request Tile::recv_data( mpi::communicator& comm, int orig, int tag)
{
  //std::cout << "RECV from " << orig << "\n";
  Mesh& mesh = get_data(); 

  mpi::request req;
  req = comm.irecv(orig, cid, mesh.mesh.data(), mesh.size() );

  return req;
}










void Solver::solve(Tile& tile) {
  Mesh& m    = tile.get_data();
  Mesh& mnew = tile.get_new_data();
  mnew.clear();

  for(int i=0; i<(int)m.Nx; i++) { 
    for(int j=0; j<(int)m.Ny; j++) { 

      // count how many is alive
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
      } else {
        mnew(i,j) = m(i,j);
      }

      // } else if(alive != 2) {
      //   mnew(i,j) = 0;
    }
  }
}



