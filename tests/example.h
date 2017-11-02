#pragma once

#include "../cell.h"


namespace example {

class Welsh : public corgi::Cell {

  public:
    Welsh(size_t i, size_t j, 
             int o, 
             size_t nx, size_t ny
             ) : corgi::Cell(i, j, o, nx, ny) { }

    // extend the base class
    std::string bark();

};


class Pembroke : public corgi::Cell {

  public:
    Pembroke(size_t i, size_t j, 
             int o, 
             size_t nx, size_t ny
             ) : corgi::Cell(i, j, o, nx, ny) { }

    // extend the base class
    std::string bark();

    // specialize the this class
    std::string howl();

};




} // end of namespace example
