#include <string>
#include <array>

#include "prtcls.h"
#include "container.h"


using namespace prtcls;


void Tile::check_outgoing_particles()
{
  for(auto&& container : containers) {
    container.check_outgoing_particles(mins, maxs);
  }
}

void Tile::delete_transferred_particles()
{
  for(auto&& container : containers) {
    container.delete_transferred_particles();
  }
}

void Tile::get_incoming_particles(
    corgi::Node<2>& grid)
{

  std::array<double,3> grid_mins = {
    static_cast<double>( grid.get_xmin() ),
    static_cast<double>( grid.get_ymin() ),
    static_cast<double>( 0.0             )
  };

  std::array<double,3> grid_maxs = {
    static_cast<double>( grid.get_xmax() ),
    static_cast<double>( grid.get_ymax() ),
    static_cast<double>( 1.0             )
  };

  // fetch incoming particles from neighbors around me
  int k = 0;
  for(int i=-1; i<=1; i++) {
    for(int j=-1; j<=1; j++) {

      // get neighboring tile
      auto ind = this->neighs(i, j); 
      uint64_t cid = 
      grid.id( std::get<0>(ind), std::get<1>(ind) );
      Tile& external_tile = 
        dynamic_cast<Tile&>( grid.get_tile(cid) );

      // loop over all containers
      for(size_t ispc=0; ispc<Nspecies(); ispc++) {
        ParticleBlock& container = get_container(ispc);
        ParticleBlock& neigh = external_tile.get_container(ispc);

        container.transfer_and_wrap_particles(
            neigh, {i,j,k}, grid_mins, grid_maxs);
      }
    }
  }

  return;
}


// create MPI tag given tile id and extra layer of differentiation
int get_tag(int cid, int extra_param)
{
  assert(extra_param < 100);
  return cid + extra_param*1e6;
}

// create MPI tag for extra communication given tile 
// id and extra layer of differentiation
int get_extra_tag(int cid, int extra_param)
{
  assert(extra_param < 100);
  return cid + extra_param*1e8;
}


void Tile::send_data( 
    mpi4cpp::mpi::communicator& comm, int dest, int tag,
    std::vector<mpi4cpp::mpi::request>& reqs
    )
{
  if(tag == 0)      return Tile::send_particle_data(comm, dest, reqs);
  else if(tag == 1) return Tile::send_particle_extra_data(comm,dest,reqs);
  else assert(false);
}

void
Tile::send_particle_data( 
    mpi4cpp::mpi::communicator& comm, int dest,
    std::vector<mpi4cpp::mpi::request>& reqs 
    )
{
  for(size_t ispc=0; ispc<Nspecies(); ispc++) {
    ParticleBlock& container = get_container(ispc);

    reqs.emplace_back(
        comm.isend(dest, get_tag(cid, ispc), 
          container.outgoing_particles.data(), 
          container.outgoing_particles.size())
        );
  }
}


void
Tile::send_particle_extra_data( 
    mpi4cpp::mpi::communicator& comm, int dest,
    std::vector<mpi4cpp::mpi::request>& reqs
    )
{
  for(size_t ispc=0; ispc<Nspecies(); ispc++) {
    ParticleBlock& container = get_container(ispc);

    if(!container.outgoing_extra_particles.empty()) {
      reqs.emplace_back(
          comm.isend(dest, get_extra_tag(cid, ispc), 
            container.outgoing_extra_particles.data(), 
            container.outgoing_extra_particles.size())
          );
    }
  }
}


void
Tile::recv_data( 
    mpi4cpp::mpi::communicator& comm, int orig, int tag,
    std::vector<mpi4cpp::mpi::request>& reqs
)
{
  if(tag == 0)      return Tile::recv_particle_data(comm,orig, reqs);
  else if(tag == 1) return Tile::recv_particle_extra_data(comm,orig, reqs);
  else assert(false);
}


void
Tile::recv_particle_data( 
    mpi4cpp::mpi::communicator& comm, int orig,
    std::vector<mpi4cpp::mpi::request>& reqs 
    )
{
  for (size_t ispc=0; ispc<Nspecies(); ispc++) {
    ParticleBlock& container = get_container(ispc);
    container.incoming_particles.resize( container.optimal_message_size );

    reqs.push_back(
        comm.irecv(orig, get_tag(cid, ispc),
          container.incoming_particles.data(),
          container.optimal_message_size)
        );
  }
}


void
Tile::recv_particle_extra_data( 
    mpi4cpp::mpi::communicator& comm, int orig,
    std::vector<mpi4cpp::mpi::request>& reqs
    )
{
  // this assumes that wait for the first message is already called
  // and passed.

  int extra_size;
  for (size_t ispc=0; ispc<Nspecies(); ispc++) {
    ParticleBlock& container = get_container(ispc);
    InfoParticle msginfo(container.incoming_particles[0]);

    // check if we need to expect extra message
    extra_size = msginfo.size() - container.optimal_message_size;
    if(extra_size > 0) {
      container.incoming_extra_particles.resize(extra_size);

      reqs.emplace_back(
          comm.irecv(orig, get_extra_tag(cid, ispc),
            container.incoming_extra_particles.data(),
            extra_size)
          );
    } else {
      container.incoming_extra_particles.clear();
    }

    //TODO: dynamic optimal_message_size here
    //container.optimal_message_size = msginfo.size();

  }
}


void Tile::pack_outgoing_particles()
{
  for(auto&& container : containers) {
    container.pack_outgoing_particles();
  }
}


void Tile::unpack_incoming_particles()
{
  for(auto&& container : containers) {
    container.unpack_incoming_particles();
  }
}


void Tile::delete_all_particles()
{
  for(auto&& container : containers) {
    container.resize(0);
  }
}



void Pusher::solve(Tile& tile) 
{
  for(auto&& container : tile.containers) {
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





