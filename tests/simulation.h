#pragma once

#include "../corgi.h"


namespace example {

class Welsh : public corgi::Cell {

  public:
    Welsh(size_t i, size_t j, int o) : corgi::Cell(i, j, o) { }

    // extend the base class
    std::string bark();

};


class Pembroke : public corgi::Cell {

  public:
    Pembroke(size_t i, size_t j, int o) : corgi::Cell(i, j, o) { }

    // extend the base class
    std::string bark();

};




} // end of namespace example
