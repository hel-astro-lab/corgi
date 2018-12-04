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


}



