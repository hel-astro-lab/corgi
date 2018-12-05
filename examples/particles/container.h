#pragma once

#include <vector>
#include <map>


namespace prtcls {


class ParticleBlock {

  protected:

  size_t Nprtcls;

  std::vector< std::vector<double> > locArr;
  std::vector< std::vector<double> > velArr;
  std::vector<double> wgtArr;

  public:

  //! multimap of particles going to other tiles
  typedef std::multimap<std::tuple<int,int,int>, int> mapType;
  mapType to_other_tiles;

  // size of the internal mesh
  size_t Nx;
  size_t Ny;
  size_t Nz;

  /// Constructor 
  ParticleBlock(size_t Nx, size_t Ny, size_t Nz);
    

  // default virtual dtor
  virtual ~ParticleBlock() = default;

  //--------------------------------------------------
    
  /// reserve memory for particles
  virtual void reserve(size_t N);

  // resize everything
  virtual void resize(size_t N);

  /// size of the container (in terms of particles)
  size_t size();

  //--------------------------------------------------
  // locations
  virtual inline double loc( size_t idim, size_t iprtcl ) const
  {
    return locArr[idim][iprtcl];
  }

  virtual inline double& loc( size_t idim, size_t iprtcl )       
  {
    return locArr[idim][iprtcl];
  }

  virtual inline std::vector<double>& loc(size_t idim) 
  {
    return locArr[idim];
  }

  //--------------------------------------------------
  // velocities
  virtual inline double vel( size_t idim, size_t iprtcl ) const
  {
    return velArr[idim][iprtcl];
  }

  virtual inline double& vel( size_t idim, size_t iprtcl )       
  {
    return velArr[idim][iprtcl];
  }

  virtual inline std::vector<double>& vel(size_t idim) 
  {
    return velArr[idim];
  }

  //--------------------------------------------------
  // weights
  virtual inline double wgt( size_t iprtcl ) const
  {
    return wgtArr[iprtcl];
  }

  virtual inline double& wgt( size_t iprtcl )       
  {
    return wgtArr[iprtcl];
  }

  virtual inline std::vector<double>& wgt() 
  {
    return wgtArr;
  }

  // --------------------------------------------------
  // particle creation & destruction methods

  virtual void add_particle (
      std::vector<double> prtcl_loc,
      std::vector<double> prtcl_vel,
      double prtcl_wgt);

  // --------------------------------------------------

  /// check and mark particles exceeding given limits
  void check_outgoing_particles(
      std::array<double,2>&,
      std::array<double,2>& );


  /// delete particles that went beyond boundaries, i.e.,
  // ended up in to_other_tiles box
  void delete_transferred_particles();


  /// transfer particles between blocks
  void transfer_and_wrap_particles(
      ParticleBlock&, 
      std::array<int,3>,
      std::array<double,3>&,
      std::array<double,3>&);

};




}
