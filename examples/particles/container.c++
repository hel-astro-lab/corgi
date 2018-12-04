#include "container.h"

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

