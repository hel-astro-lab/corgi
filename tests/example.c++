#include <string>

#include "example.h"

using namespace example;


std::string Welsh::bark() { return "Woof!"; };

std::string Pembroke::bark() { return "Ruff!"; };
std::string Pembroke::howl() { return "Auuuuuu!"; };


