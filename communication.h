#pragma once

#include <mpi4cpp/mpi.h>


// introduce class to mpi
namespace mpi4cpp { namespace mpi {

template <>
struct is_mpi_datatype<corgi::Communication> 
  : mpl::true_ { };

template <>
inline MPI_Datatype get_mpi_datatype<corgi::Communication>(
    const corgi::Communication& obj) 
{
  MPI_Datatype obj_type;
    
  //--------------------------------------------------
  // check if we have already created the object once
  std::type_info const* t = &typeid(::corgi::Communication);
  obj_type = mpi4cpp::mpi::detail::mpi_datatype_cache().get(t);

  if(obj_type != MPI_DATATYPE_NULL) return obj_type;

  //--------------------------------------------------
  // if not; we'll create it now

  // get start of the class
  MPI_Aint base;
  MPI_Get_address( &obj, &base ); 

  //--------------------------------------------------
  // Communication class interiors:
  //
  // 1 int cid;
  // 2 std::array<int, 3> indices;
  // 3 int owner;
  // 4 int top_virtual_owner;
  // 5 int communications = 0;
  // 6 int number_of_virtual_neighbors = 0;
  // 7 std::array<double, 3> mins;
  // 8 std::array<double, 3> maxs;
  // x bool local = false;

  // how many elements per each type
  std::array<int, 8> block_lengths{
    { 1, 3, 1, 1, 1, 1, 3, 3}
  };

  // and then the actual members
  std::array<MPI_Aint, 8> member_offsets; // relative offsets
  MPI_Get_address( &obj.cid,                         &member_offsets[0]);
  MPI_Get_address( &obj.indices[0],                  &member_offsets[1]);
  MPI_Get_address( &obj.owner,                       &member_offsets[2]);
  MPI_Get_address( &obj.top_virtual_owner,           &member_offsets[3]);
  MPI_Get_address( &obj.communications,              &member_offsets[4]);
  MPI_Get_address( &obj.number_of_virtual_neighbors, &member_offsets[5]);
  MPI_Get_address( &obj.mins[0],                     &member_offsets[6]);
  MPI_Get_address( &obj.maxs[0],                     &member_offsets[7]);

  // create real (absolute) offsets (=rel - base)
  std::array<MPI_Aint, 8> offsets {
    member_offsets[0] - base,
    member_offsets[1] - base,
    member_offsets[2] - base,
    member_offsets[3] - base,
    member_offsets[4] - base,
    member_offsets[5] - base,
    member_offsets[6] - base,
    member_offsets[7] - base
  };

  /*
  std::cout << "offsets are "
    << " 0 " << offsets[0] << " / " << member_offsets[0] << " x " << base 
    << " 1 " << offsets[1] << " / " << member_offsets[1] << " x " << base 
    << " 2 " << offsets[2] << " / " << member_offsets[2] << " x " << base 
    << "\n";
  */

  // introduce datatypes
  std::array<MPI_Datatype, 8> datatypes{
    { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE }
  };

  //--------------------------------------------------
  // create datatype; this is standard format and should not be changed
  MPI_Type_create_struct(
      block_lengths.size(),
      block_lengths.data(),
      offsets.data(),
      datatypes.data(),
      &obj_type);

  MPI_Type_commit(&obj_type);

  // check and correct for extent in arrays
  MPI_Aint lb, extent, dist[2];
  MPI_Type_get_extent(obj_type, &lb, &extent);
  std::array<corgi::Communication,10> arr;

  MPI_Get_address(&arr[0], &dist[0]);
  MPI_Get_address(&arr[1], &dist[1]);

  MPI_Datatype temptype;
  if (extent != (dist[1] - dist[0])) {
    std::cout << "WARNING! Extra padding between arrays in MPI datatype\n";
    temptype = obj_type;
    lb = 0;
    extent = dist[1] - dist[0];

    MPI_Type_create_resized(temptype, lb, extent, &obj_type);
    MPI_Type_commit(&obj_type);
    MPI_Type_free(&temptype);
  }
  
  // save datatype to internal cache
  mpi4cpp::mpi::detail::mpi_datatype_cache().set(t, obj_type);

  return obj_type;
}



} } // ns mpi4cpp::mpi
