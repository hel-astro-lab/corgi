#include <string>

#include "example.h"

using namespace example;


// Welsh methods
std::string Welsh::bark() { return "Woof!"; };

// Pembroke methods
std::string Pembroke::bark() { return "Ruff!"; };
std::string Pembroke::howl() { return "Auuuuuu!"; };


// Vallhund definitions
std::string Viking::prayForOdin() { return "Oooh Mighty Odin!"; };
std::string Swede::fika() { return "---: It is fika time, get the kanelbullas"; };
std::string Vallhund::bark() { return "ruf ruf ruf"; 
};

// Grid methods
std::string Grid::petShop() { return "No Corgis for sale."; };


