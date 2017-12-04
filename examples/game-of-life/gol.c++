#include <string>

#include "gol.h"
#include "../../toolbox/dataContainer.h"


using namespace gol;

/// initialize internal mesh
Mesh::Mesh(int Nx, int Ny) : Nx(Nx), Ny(Ny) {
  mesh.resize((Nx + 2*halo) * (Ny + 2*halo));    
}


/// Copy vertical slice
void Mesh::copyVert(Mesh& rhs, int lhsI, int rhsI) {
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int j=0; j<this->Ny; j++) { 
    this->operator()(lhsI, j) = rhs(rhsI, j);
  }
}


/// Copy horizontal slice 
void Mesh::copyHorz(Mesh& rhs, int lhsJ, int rhsJ) {
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");

  for(int i=0; i<this->Nx; i++) { 
    this->operator()(i, lhsJ) = rhs(i, rhsJ);
  }
}


// --------------------------------------------------


/// Add data to the container
void CellularAutomataCell::addData(Mesh m) {
  data.push_back(m);
}


/// get current patch
Mesh& CellularAutomataCell::getData() {
  return *data.get();
};

Mesh* CellularAutomataCell::getDataPtr() {
  return data.get();
};

/// get new data
Mesh& CellularAutomataCell::getNewData() {
  return *data.getNew();
};


/// Update boundary/halo regions from neighbors
void CellularAutomataCell::updateBoundaries(Grid& grid) {

  Mesh& mesh = getData(); // target as a reference to update into
  

  // left 
  CellPtr cleft = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( neighs(-1, 0) ));
  Mesh& mleft = cleft->getData();
  mesh.copyVert(mleft, -1, mleft.Nx-1); // copy from right side to left

  // right
  CellPtr cright = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( neighs(+1, 0) ));
  Mesh& mright = cright->getData();
  mesh.copyVert(mright, mesh.Nx, 0); // copy from left side to right

  // top 
  CellPtr ctop = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( neighs(0, +1) ));
  Mesh& mtop = ctop->getData();
  mesh.copyHorz(mtop, mesh.Ny, 0); // copy from bottom side to top

  // bottom
  CellPtr cbot = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( neighs(0, -1) ));
  Mesh& mbot = cbot->getData();
  mesh.copyHorz(mbot, -1, mbot.Ny-1); // copy from top side to bottom




  // diagonals/corners
  // --------------------------------------------------  
    
  // top right
  CellPtr ctr = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( this->neighs(+1, +1) ));
  Mesh& mtr = ctr->getData();
  mesh(mesh.Nx, mesh.Ny) = mtr(0,0);
  

  // top left
  CellPtr ctl = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( this->neighs(-1, +1) ));
  Mesh& mtl = ctl->getData();
  mesh(-1, mesh.Ny) = mtl(mtl.Nx-1,0);


  // bottom right
  CellPtr cbr = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( this->neighs(+1, -1) ));
  Mesh& mbr = cbr->getData();
  mesh(mesh.Nx, -1) = mbr(0,mbr.Ny-1);

  // bottom left WORKS
  CellPtr cbl = std::dynamic_pointer_cast<CellType>(grid.getCellPtr( this->neighs(-1, -1) ));
  Mesh& mbl = cbl->getData();
  mesh(-1,-1) = mbl(mbl.Nx-1, mbl.Ny-1);
  

};


void Solver::solve(CellularAutomataCell& cell) {
  Mesh& m    = cell.getData();
  Mesh& mnew = cell.getNewData();


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



// std::string Grid::petShop() { return "No Corgis for sale."; };

