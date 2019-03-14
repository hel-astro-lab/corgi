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
void Tile::add_data(Mesh m) {
  data.push_back(m);
}


/// get current patch
Mesh& Tile::get_data() {
  return data.get();
}

/// get new data
Mesh& Tile::get_new_data() {
  return data.get(1);
}


/// Update boundary/halo regions from neighbors
void Tile::update_boundaries(corgi::Node<2>& grid) 
{
  int ito, jto, ifro, jfro;
  Tileptr tpr;

  Mesh& mesh = get_data(); // target as a reference to update into

  for(int in=-1; in <= 1; in++) {
    for(int jn=-1; jn <= 1; jn++) {
      if (in == 0 && jn == 0) continue;

      tpr = std::dynamic_pointer_cast<Tile_t>(grid.get_tileptr( neighs(in, jn) ));
      if (tpr) {
        Mesh& mpr = tpr->get_data();

        /* diagonal rules are:
        if + then to   n
        if + then from 0

        if - then to   -1
        if - then from n-1
        */

        if (in == +1) { ito = mesh.Nx; ifro = 0; }
        if (jn == +1) { jto = mesh.Ny; jfro = 0; }

        if (in == -1) { ito = -1;      ifro = mpr.Nx-1; }
        if (jn == -1) { jto = -1;      jfro = mpr.Ny-1; }

        // copy
        if      (jn == 0) mesh.copy_vert(mpr, ito, ifro);   // vertical
        else if (in == 0) mesh.copy_horz(mpr, jto, jfro);   // horizontal
        else              mesh(ito, jto) = mpr(ifro, jfro); // diagonal
        
      } // end of if(tpr)
    }
  }

}

std::vector<mpi4cpp::mpi::request> Tile::send_data( mpi4cpp::mpi::communicator& comm, int dest, int /*tag*/)
{
  //std::cout << "SEND to " << dest << "\n";
  Mesh& mesh = get_data(); 

  std::vector<mpi4cpp::mpi::request> reqs;
  reqs.push_back( comm.isend(dest, cid, mesh.mesh.data(), mesh.size()) );

  return reqs;
}

std::vector<mpi4cpp::mpi::request> Tile::recv_data( mpi4cpp::mpi::communicator& comm, int orig, int /*tag*/)
{
  //std::cout << "RECV from " << orig << "\n";
  Mesh& mesh = get_data(); 

  std::vector<mpi4cpp::mpi::request> reqs;
  reqs.push_back( comm.irecv(orig, cid, mesh.mesh.data(), mesh.size()) );

  return reqs;
}


void Solver::solve(Tile& tile) {
  Mesh& m    = tile.get_data();
  Mesh& mnew = tile.get_new_data();
  mnew.clear();

  for(int i=0; i<(int)m.Nx; i++) { 
    for(int j=0; j<(int)m.Ny; j++) { 

      // count how many are alive
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



