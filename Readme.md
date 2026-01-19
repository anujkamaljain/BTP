# Opportunistic Routing Simulation

This project implements and simulates opportunistic routing algorithms in C++.

## Phases

### Phase 1 (Phase1.cpp)
Implements the core opportunistic routing algorithm. It computes routing tables and expected transmission costs for reaching sets of destinations using forwarding sets.

### Phase 2 (phase2.cpp)
Extends Phase 1 by adding packet transmission simulation. It simulates multiple packet transmissions to validate theoretical costs with empirical results, showing average, min, and max transmissions.

### Phase 3 (phase3.cpp)
Further extends Phase 2 by incorporating network coding simulation. It demonstrates the benefits of network coding, such as XOR transmissions, for efficient multicast to multiple destinations.

## Usage
Compile and run the respective .cpp files with appropriate inputs for network topology, source, destinations, and simulation parameters.

“All rights reserved. Unauthorized use prohibited.”