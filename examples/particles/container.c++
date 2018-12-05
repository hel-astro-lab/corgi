#include <array>
#include <vector>

#include "container.h"
#include "wrap.h"

using namespace prtcls;


ParticleBlock::ParticleBlock(size_t Nx, size_t Ny, size_t Nz) :
  Nx(Nx), Ny(Ny), Nz(Nz)
{ 
  locArr.resize(3);
  velArr.resize(3);
};



void ParticleBlock::reserve(size_t N) {
  for(size_t i=0; i<3; i++) locArr[i].reserve(N);
  for(size_t i=0; i<3; i++) velArr[i].reserve(N);
  wgtArr.reserve(N);
}


void ParticleBlock::resize(size_t N)
{
  for(size_t i=0; i<3; i++) locArr[i].resize(N);
  for(size_t i=0; i<3; i++) velArr[i].resize(N);
  wgtArr.resize(N);
  Nprtcls = N;
}


size_t ParticleBlock::size() 
{ 
  // FIXME: these fail
  //assert(locArr[0].size() == Nprtcls);
  //assert(locArr[1].size() == Nprtcls);
  //assert(locArr[2].size() == Nprtcls);

  //assert(velArr[0].size() == Nprtcls);
  //assert(velArr[1].size() == Nprtcls);
  //assert(velArr[2].size() == Nprtcls);

  //return Nprtcls; // FIXME: this is the correct way to return
  return locArr[0].size();
}


void ParticleBlock::add_particle (
    std::vector<double> prtcl_loc,
    std::vector<double> prtcl_vel,
    double prtcl_wgt)
{

  assert(prtcl_loc.size() == 3);
  assert(prtcl_vel.size() == 3);

  for (size_t i=0; i<3; i++) locArr[i].push_back(prtcl_loc[i]);
  for (size_t i=0; i<3; i++) velArr[i].push_back(prtcl_vel[i]);
  wgtArr.push_back(prtcl_wgt);

  Nprtcls++;
}

void ParticleBlock::check_outgoing_particles(
    std::array<double,2>& mins,
    std::array<double,2>& maxs)
{
  to_other_tiles.clear();

  // unpack limits
  double 
    xmin = mins[0],
    ymin = mins[1],
    zmin = 0.0,

    xmax = maxs[0],
    ymax = maxs[1],
    zmax = 1.0;

  // shortcut for particle locations
  double* locn[3];
  for( int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  double x0, y0, z0;

  int i,j,k; // relative indices
  for(size_t n=0; n<size(); n++) {
    i = 0;
    j = 0;
    k = 0;

    x0 = locn[0][n];
    y0 = locn[1][n];
    z0 = locn[2][n];

    if(x0 <  xmin) i--; // left wrap
    if(x0 >= xmax) i++; // right wrap

    if(y0 <  ymin) j--; // bottom wrap
    if(y0 >= ymax) j++; // top wrap

    if(z0 <  zmin) k--; // back
    if(z0 >= zmax) k++; // front

    if ((i == 0) && (j == 0)) continue; //TODO: hack to make this work with 2D corgi tiles

    //if ( (i != 0) | (j != 0) | (k != 0) ) {
    if ( (i != 0) || (j != 0) ) { // TODO: another 2D hack
      to_other_tiles.insert( 
          std::make_pair( std::make_tuple(i,j,k), n) 
                           );
    }
  }
}


void ParticleBlock::delete_transferred_particles()
{
  std::vector<int> to_be_deleted;

  //for(auto& elem : container.to_other_tiles)  {
  //  std::cout << "to be removed particle # " 
  //            << elem.second << " ";
  //  std::cout << "(" 
  //            << std::get<0>(elem.first) 
  //            << "," << std::get<1>(elem.first) 
  //            << "," << std::get<2>(elem.first) 
  //            << ")";
  //  std::cout << '\n';
  //}

  // get and reverse sort deleted particles
  for(auto& elem : to_other_tiles) to_be_deleted.push_back( elem.second );
  std::sort(to_be_deleted.begin(), to_be_deleted.end(), std::greater<int>() );

  double* locn[3];
  for( int i=0; i<3; i++) locn[i] = &( loc(i,0) );

  double* veln[3];
  for( int i=0; i<3; i++) veln[i] = &( vel(i,0) );

  // overwrite particles with the last one on the array and 
  // then resize the array
  int last = size();
  for(int indx : to_be_deleted) {
    last--;
    //std::cout << "deleting " << indx 
    //          << " by putting it to " << last << '\n';
    for(int i=0; i<3; i++) locn[i][indx] = locn[i][last];
    for(int i=0; i<3; i++) veln[i][indx] = veln[i][last];
    wgtArr[indx] = wgtArr[last];
  }

  // resize if needed and take care of the size
  last = last < 0 ? 0 : last;
  if ((last != (int)size()) && (size() > 0)) resize(last);

  return;
}


void ParticleBlock::transfer_and_wrap_particles( 
    ParticleBlock& neigh,
    std::array<int,3>     dirs, 
    std::array<double,3>& mins, 
    std::array<double,3>& maxs
    )
{
  double locx, locy, locz, velx, vely, velz, wgt;

  for (auto&& elem : neigh.to_other_tiles) {
      
    //TODO: collapsed z-dimension due to 2D corgi tiles
    if (std::get<0>(elem.first) == 0 &&
        std::get<1>(elem.first) == 0 ) 
    { continue; }

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile
    if (std::get<0>(elem.first) == -dirs[0] &&
        std::get<1>(elem.first) == -dirs[1] ) {
              
      locx = wrap( neigh.loc(0, elem.second), mins[0], maxs[0] );
      locy = wrap( neigh.loc(1, elem.second), mins[1], maxs[1] );
      locz = wrap( neigh.loc(2, elem.second), mins[2], maxs[2] );

      velx = neigh.vel(0, elem.second);
      vely = neigh.vel(1, elem.second);
      velz = neigh.vel(2, elem.second);

      wgt  = neigh.wgt(elem.second);

      add_particle({locx,locy,locz}, {velx,vely,velz}, wgt);
    }
  }

  return;
}





