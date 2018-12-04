#include <string>

#include "prtcls.h"
#include "container.h"


using namespace prtcls;
using namespace mpi4cpp;




mpi::request Tile::send_data( mpi::communicator& comm, int dest, int tag)
{
  //std::cout << "SEND to " << dest << "\n";

  mpi::request req;
  //req = comm.isend(dest, cid, mesh.mesh.data(), mesh.size() );

  return req;
}

mpi::request Tile::recv_data( mpi::communicator& comm, int orig, int tag)
{
  //std::cout << "RECV from " << orig << "\n";

  mpi::request req;
  //req = comm.irecv(orig, cid, mesh.mesh.data(), mesh.size() );

  return req;
}


void Pusher::solve(Tile& tile) 
{

  for (size_t ispc=0; ispc<tile.Nspecies(); ispc++) {
    ParticleBlock& container = tile.get_container(ispc);

    int nparts = container.size();

    // initialize pointers to particle arrays
    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( container.vel(i,0) );


    // loop over particles
    int n1 = 0;
    int n2 = nparts;

    double u0, v0, w0;
    double g;
    double c = 0.5;

      for(int n=n1; n<n2; n++) {

        u0 = c*vel[0][n];
        v0 = c*vel[1][n];
        w0 = c*vel[2][n];

        // position advance
        g = c / sqrt(c*c + u0*u0 + v0*v0 + w0*w0);
        for(size_t i=0; i<3; i++)
          loc[i][n] += vel[i][n]*g*c;

      }

  }//end of loop over species

}





