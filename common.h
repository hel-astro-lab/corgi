#pragma once

#include <cstdint>
#include <climits>

/*
#include <limits>
#include <vector>

#define sqr(x) ((x)*(x))
#define pow2(x) sqr(x)
#define pow3(x) ((x)*(x)*(x))
*/


#ifndef MASTER_RANK
#define MASTER_RANK 0
#endif


// Define uint64_t for MPI
//#ifndef MPI_UINT64_T
//#define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
//#endif

// Define size_t for MPI
#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "Error while initializing size_t type for MPI"
#endif


// TODO remove double definition in Node (done for python bindings)
// struct globalflags {
//     static bool master;
// };

/*
namespace conf {

    /// Grid dimensions
    size_t Nx = 10;
    size_t Ny = 10;

    /// block size inside spatial tile
    size_t NxTile = 2;
    size_t NyTile = 2;

    /// physical grid dimensions
    double xmin = 0.0;
    double xmax = 1.0;

    double ymin = 0.0;
    double ymax = 1.0;

}
*/

/*
namespace BC {

    /// Periodic x boundary condition
    size_t xwrap( int i ) {
        while (i < 0) {
            i += conf::Nx;
        }
        while (i >= conf::Nx) {
            i -= conf::Nx;
        }
        return size_t(i);
    }

    /// Periodic y boundary condition
    size_t ywrap( int j ) {
        while (j < 0) {
            j += conf::Ny;
        }
        while (j >= conf::Ny) {
            j -= conf::Ny;
        }
        return size_t(j);
    }

}
*/


/*! Contains different tile boundary types. Type of computation will depend on 
 * what this type is.
namespace tileType {
    enum {
        LOCAL,    //! Default type indicating that tile is owned by the current process
        VIRTUAL,  //! virtual tile (owned by other process)
        OUTFLOW,  //! outflow tile
        N_CELLTYPES
    };
}
*/

namespace commType {
    enum {
        NTILES,   //! Number of incoming tiles,
        TILEDATA, //! Tile data array
        N_COMMTYPES
    };
}




