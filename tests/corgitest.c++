#include <string>

#include "corgitest.h"

using namespace corgitest;


// Welsh methods
std::string Welsh::bark() { return "Woof!"; }

// Pembroke methods
std::string Pembroke::bark() { return "Ruff!"; }
std::string Pembroke::howl() { return "Auuuuuu!"; }


// Vallhund definitions
std::string Viking::pray_for_odin() { return "Oooh Mighty Odin!"; }
std::string Swede::fika() { return "---: It is fika time, get the kanelbullas"; }
std::string Vallhund::bark() { return "ruf ruf ruf"; }

// Grid methods
//std::string Grid::pet_shop() { return "No Corgis for sale."; }
