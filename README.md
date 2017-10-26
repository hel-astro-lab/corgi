# CORGI 
#### C++ Object-oRiented Grid Infastructure
##### (...or some other combination of relevant words that spell out another funny word.)
---

`corgi` is a c++ template library for building massively parallel grids. It is aims to be *completely* local, so every node will only need to work on the data that it currently sees/owns. Underlying communication is done using mpi. Novelty is in the load balancing scheme that applies cellular automata kind of rules to even out load imbalances and tries, at the same time, to minimize the inter-node communication.

**This is still work-in-progress**


## In a nutshell:
- c++-14 compatible compiler is needed,
- `mpi` library should be installed, and
- it is header only; just include and start using!



