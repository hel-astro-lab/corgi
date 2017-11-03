#pragma once

#include "../cell.h"
#include "../corgi.h"


namespace example {

class Welsh : public corgi::Cell {

  public:
    Welsh(size_t i, size_t j, 
             int o, 
             size_t nx, size_t ny
             ) : corgi::Cell(i, j, o, nx, ny) { }

    ~Welsh() { };

    // extend the base class
    std::string bark();

};


class Pembroke : public corgi::Cell {

  public:
    Pembroke(size_t i, size_t j, 
             int o, 
             size_t nx, size_t ny
             ) : corgi::Cell(i, j, o, nx, ny) { }

    ~Pembroke() { };

    // extend the base class
    std::string bark();

    // specialize the this class
    std::string howl();

};


class Grid : public corgi::Node {

  public:
    Grid(size_t nx, size_t ny) : corgi::Node(nx, ny) { }

    ~Grid() { };

    // std::unordered_map< uint64_t, std::shared_ptr<corgi::Cell>> cells;

    std::string petShop();


};



} // end of namespace example
