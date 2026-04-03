# Opportunistic Routing with Network Coding — Simulation Framework

A progressive C++ simulation of **opportunistic routing** and **network coding** techniques for wireless networks. The project is organized into multiple phases (Phase 1 through Phase 7a), each building upon the previous one to add more realistic and complex capabilities.

---

## Table of Contents

- [Overview](#overview)
- [Core Concepts](#core-concepts)
- [Phase Summary](#phase-summary)
- [Phase 1 — Core Routing Algorithm](#phase-1--core-routing-algorithm)
- [Phase 2 — Packet Transmission Simulation](#phase-2--packet-transmission-simulation)
- [Phase 3 — Network Coding (XOR)](#phase-3--network-coding-xor)
- [Phase 4 — Multi-Source NC with Coding Window](#phase-4--multi-source-nc-with-coding-window)
- [Phase 5 — Multi-Source Multi-Destination (2-Hop)](#phase-5--multi-source-multi-destination-2-hop)
- [Phase 6 — Multi-Source Multi-Destination Multi-Hop (Optimal FS)](#phase-6--multi-source-multi-destination-multi-hop-optimal-fs)
- [Phase 6a — Multi-Hop with Heuristic (Singleton FS)](#phase-6a--multi-hop-with-heuristic-singleton-fs)
- [Phase 7 — Optimized Multi-Hop (Optimal FS + Performance)](#phase-7--optimized-multi-hop-optimal-fs--performance)
- [Phase 7a — Optimized Multi-Hop (Heuristic Singleton FS + Performance)](#phase-7a--optimized-multi-hop-heuristic-singleton-fs--performance)
- [How to Build & Run](#how-to-build--run)
- [Input File Format](#input-file-format)
- [Project Structure](#project-structure)

---

## Overview

This project simulates **opportunistic routing** in wireless networks — a technique where a packet's next hop is decided **after** transmission, based on which neighbors actually received it. This is combined with **network coding** (XOR-based), where coded packets allow multiple destinations to recover missing data simultaneously, reducing overall transmissions.

## Core Concepts

| Concept | Description |
|---|---|
| **ETX (Expected Transmission Count)** | Edge cost = avg transmissions needed for one successful delivery. Probability of success = `1 / cost`. |
| **Forwarding Set (FS)** | Subset of neighbors that a node transmits to. The algorithm finds the optimal FS that minimizes expected cost. |
| **Routing Table** | Stores the minimum expected cost from each node to every destination subset, along with the best forwarding set. |
| **Iterative Cost Calculation** | Uses fixed-point iteration over all `2^n` reception outcomes (bitmask enumeration) to compute expected transmission cost. Tolerance `ε = 0.001`, max 1000 iterations. |
| **Coding Window** | A buffer at forwarding nodes. Packets fill the window; coded packets are transmitted. Once all destinations receive `≥ window size` coded packets, all original packets are decoded/delivered. The window is then cleared and refilled. |
| **XOR Recovery** | After Phase 1 delivery, missing packets at destinations are recovered via pairwise XOR or individual retransmission. |
| **L'Ecuyer RNG** | Phases 4–7a use the L'Ecuyer combined multiplicative congruential generator (`randomgeneration()`) for reproducible, high-quality random numbers, seeded with `-(long)time(NULL)`. |
| **TeeBuf** | Phases 5–7a implement a `TeeBuf` class (custom `streambuf`) that duplicates all `cout` output to both the console and an output `.txt` file simultaneously. |

---

## Phase Summary

| Phase | File | Lines | Input | Graph | Sources | Destinations | Hops | Network Coding | Coding Window | RNG | Forwarding Set Selection | Output Tee |
|-------|------|-------|-------|-------|---------|--------------|------|----------------|---------------|-----|--------------------------|------------|
| 1 | `Phase1.cpp` | 395 | stdin | Directed | 1 | 1–2 | Any | ✗ | ✗ | — | Singletons + full set | ✗ |
| 2 | `phase2.cpp` | 652 | `phase2_input.txt` | Directed | 1 | 1–2 | Any | ✗ | ✗ | `mt19937` | Singletons + full set | ✗ |
| 3 | `phase3.cpp` | 736 | stdin | Directed | 1 | 2 (required) | Any | ✓ (XOR) | ✗ | `mt19937` | Singletons + full set | ✗ |
| 4 | `phase4.cpp` | 1255 | `phase4_input.txt` | Directed | 1–2 | 2 | 2 | ✓ (XOR) | ✓ | L'Ecuyer | Singletons + full set | ✗ |
| 5 | `phase5.cpp` | 1035 | `phase5_input.txt` | **Undirected** | 2–4 | 2+ each | 2 | ✓ (XOR) | ✓ | L'Ecuyer | **All subsets** | ✓ |
| 6 | `phase6.cpp` | 1134 | `phase6_input.txt` | **Undirected** | 3–6 | 3 each | **3+** | ✓ (XOR) | ✓ | L'Ecuyer | **All subsets** | ✓ |
| 6a | `phase6a.cpp` | 1149 | `phase6a_input.txt` | **Undirected** | 3–6 | 3 each | **3+** | ✓ (XOR) | ✓ | L'Ecuyer | **Heuristic: singleton FS for simulation** | ✓ |
| 7 | `phase7.cpp` | 1132 | `phase7_input.txt` | **Undirected** | 3–6 | 3 each | **3+** | ✓ (XOR) | ✓ | L'Ecuyer | **All subsets** | ✓ |
| 7a | `phase7a.cpp` | 1147 | `phase7a_input.txt` | **Undirected** | 3–6 | 3 each | **3+** | ✓ (XOR) | ✓ | L'Ecuyer | **Heuristic: singleton FS for simulation** | ✓ |

---

## Phase 1 — Core Routing Algorithm

**File:** `Phase1.cpp` (395 lines) · **Input:** stdin (interactive prompts)

### What It Does

Implements the foundational opportunistic routing algorithm. Given a directed network graph, a source, and a set of destinations, it computes:

1. **Routing tables** — expected transmission cost from every relevant node to every destination subset.
2. **Best forwarding sets** — the optimal set of neighbors to forward to for minimum cost.

### Key Code Details

- **Class:** `OpportunisticRouting` with a single `source` and `set<int> destinations`.
- **Data Structures:**
  - `graph`: `map<int, map<int, double>>` — adjacency list with edge costs.
  - `routingTable`: `map<int, map<set<int>, double>>` — node → destination-subset → expected cost.
  - `bestForwardingSet`: `map<int, map<set<int>, set<int>>>` — node → destination-subset → best FS.
- **Cost Computation** (`calculateCostIterative`): For each candidate forwarding set, enumerates all `2^n` possible reception outcomes (success/failure per neighbor) using bitmask iteration. Computes the probability-weighted expected transmissions and converges via fixed-point iteration (initial guess `x = 5.0`, tolerance `ε = 0.001`, max 1000 iterations). Reception probability is computed as `p = 1 / cost`.
- **Forwarding Set Candidates:** Only singletons (each neighbor individually) and the full neighbor set are evaluated. This is a limitation corrected in Phase 5+.
- **Routing Table Building** (`solve`): Three-step bottom-up approach:
  1. Initialize destination nodes with cost = 0.
  2. Build routing tables for predecessors of destinations.
  3. Build routing table for the source.
- **Output:** Prints network topology, per-node routing tables with forwarding sets, and the final cost + best forwarding set from source to destinations.

### Input Format (stdin)

```
5          ← number of nodes
0          ← source node
2          ← number of destinations
3 4        ← destination node IDs
6          ← number of edges
0 1 4      ← edge: from=0, to=1, cost=4
0 2 8
1 3 5
1 4 5.5
2 3 8
2 4 8
```

---

## Phase 2 — Packet Transmission Simulation

**File:** `phase2.cpp` (652 lines) · **Input:** `phase2_input.txt`

### What's New (vs Phase 1)

- Reads input from a **text file** with `#` comment line support and a `getNextDataLine()` helper function, instead of interactive stdin.
- Adds a **Monte Carlo packet simulation** that sends thousands of packets through the network to **empirically validate** the theoretical routing costs.
- Introduces the `<random>` header and uses **`mt19937`** (Mersenne Twister) RNG with `random_device` seeding.

### Key Code Additions

- **`TransmissionOutcome` struct**: Contains `successfulNodes`, `probability`, and cumulative probability range (`cumProbStart`, `cumProbEnd`).
- **`calculateProbabilityRegions()`**: Maps the `[0, 1)` random number space to each of the `2^n` reception outcomes using cumulative probability ranges.
- **`simulatePacketTransmission()`**: Recursive function. For each packet, generates a random number, determines which neighbors received based on cumulative probability regions, removes covered destinations, and recursively continues from the best receiving node (lowest routing table cost) until all destinations are reached. Retransmits from the same node if no neighbors received.
- **`runSimulation()`**: Runs `numPackets` simulations. Reports **avg / min / max transmissions** and compares to the theoretical cost from the routing table.

### Input File Format (`phase2_input.txt`)

```
# Number of nodes
5
# Source node
0
# Number of destinations
2
# Destination nodes
3 4
# Number of edges
6
# Edges (from  to  cost)
0 1 4
0 2 8
1 3 5
1 4 5.5
2 3 8
2 4 8
# Number of packets to simulate
10000
```

### Output

Includes everything from Phase 1, plus:
- Simulation results: average, min, max transmissions
- Theoretical vs empirical cost comparison with difference

---

## Phase 3 — Network Coding (XOR)

**File:** `phase3.cpp` (736 lines) · **Input:** stdin (interactive)

### What's New (vs Phase 2)

- Introduces **network coding simulation** using XOR-based recovery.
- Adds a `simulateToAnyDestination()` method — packet is forwarded until **any one** destination receives it (not all).
- Removes verbose debug output from `addEdge()` and `addDestination()`.
- Adds `runNetworkCodingSimulation()` for batch-based NC simulation.

### Key Logic — Network Coding Simulation

The simulation operates in **two phases per batch**:

| Sub-Phase | Description |
|-----------|-------------|
| **Phase 1 (Delivery)** | Each packet is routed until **any one** destination receives it. Track which destination got which packet using `receivedPackets[dest_id] = {packet_ids}`. |
| **Phase 2 (Recovery)** | Find missing packets at each destination. Use **XOR recovery**: if dest₁ has packet A (missing at dest₂) and dest₂ has packet B (missing at dest₁), then `XOR(A, B)` satisfies both in a single coded transmission (cost = 2 hops). Remaining gaps are filled individually (also 2 hops each). |

> **Constraint:** Requires exactly **2 destinations** (hardcoded check: `if (destinations.size() != 2)`).

### Input Format (stdin)

Same as Phase 1, plus additional prompts for:
- Number of batches
- Packets per batch

### Output

- Total packets processed
- Theoretical cost to both destinations
- XOR transmissions used
- Individual transmissions used
- Average transmissions per packet **with** network coding

---

## Phase 4 — Multi-Source NC with Coding Window

**File:** `phase4.cpp` (1255 lines) · **Input:** `phase4_input.txt`

### What's New (vs Phase 3)

- **Multiple sources** (1–2) feeding packets into the network.
- **Coding window** mechanism at forwarding (intermediate) nodes.
- Replaces `std::mt19937` with the **L'Ecuyer combined multiplicative congruential RNG** (`randomgeneration()` function with macros `IM1`, `IM2`, `IA1`, etc.) for higher-quality simulation randomness. Global seed initialized as `-(long)time(NULL)`.
- File-based input with separate source and destination specifications.
- Removes the `mt19937& rng` parameter from simulation methods; uses global `global_seed` instead.
- Adds **`simulateToForwardingNode()`** and **`simulateFromForwardingToDestination()`** methods for the two-hop pipeline.
- `simulateFromForwardingToDestination()` returns `set<int>` of **all** destinations that received (not just one), enabling wireless broadcast tracking.

### Key Logic — Coding Window

```
┌─────────┐         ┌──────────────┐         ┌──────┐
│ Source 0 │──pkt──▶ │ Forwarding   │──coded──▶│Dest 4│
│ Source 1 │──pkt──▶ │ Node (2 or 3)│──coded──▶│Dest 5│
└─────────┘         └──────────────┘         └──────┘
                    ▲ Coding Window ▲
                    │  [p1,p2,...pW] │
                    └───────────────┘
```

1. **Source → Forwarding:** Sources pop packets from their queue and transmit. On success, the packet is placed into the forwarding node's queue (`forwardingQueues[node]`). On failure, retransmit.
2. **Coding Window Fill:** At each forwarding node, packets are moved from the queue into the coding window (up to `codingWindowSize`).
3. **Coded Transmission:** The forwarding node transmits a coded packet to destinations using `simulateFromForwardingToDestination()`. For each destination that receives, `windowRecvCount[node][dest]` is incremented.
4. **Window Flush:** When **every destination** has received `≥ targetWindowSize` coded packets, all packets in the window are marked as decoded at all destinations. The window is cleared and counters reset.
5. **Refill:** After flush check, the window is refilled from the queue.

### Node Selection

A random node is selected each iteration from `sources ∪ forwarding_nodes` using L'Ecuyer RNG. If it's a source, it sends a new packet. If it's a forwarding node, it processes its coding window.

### Phase 2 (XOR Recovery)

After Phase 1 completes (all sources emptied, all queues and windows flushed):
- For 2 destinations: pairwise XOR recovery with 2-hop cost each. Remaining missing packets sent individually.
- For 1 destination: individual sends only.

### Input File Format (`phase4_input.txt`)

```
# Number of nodes
6
# Number of source nodes
2
# Source nodes
0 1
# Number of destinations
2
# Destination nodes
4 5
# Number of edges
8
# Edges (from  to  cost)
0 2 4
0 3 5.5
...
# Packets per source
10000
# Coding window size
64
```

---

## Phase 5 — Multi-Source Multi-Destination (2-Hop)

**File:** `phase5.cpp` (1035 lines) · **Input:** `phase5_input.txt`

### What's New (vs Phase 4)

- **Class redesign:** Removes the single `source` and `destinations` members. Instead:
  - `srcDestNodes`: set of all nodes that act as source/destination.
  - `sourceDestMap`: `map<int, set<int>>` mapping each source to its specific destination set.
  - Constructor takes only `numNodes` (no `source` parameter).
- **Nodes act as both source and destination** simultaneously (e.g., node 0 sends to {4,5} while node 4 sends to {0,1}).
- **Undirected graph** — edges are added in both directions with the same cost in `main()`.
- Each source has its **own destination set** (source-destination mapping).
- **Forwarding queues and coding windows are maintained per `(forwarding_node, source)` pair** using `map<pair<int,int>, queue<...>>` and `map<int, map<int, vector<int>>>`, ensuring independent processing.
- **`buildRoutingTableForNode()` now accepts a `destinations` parameter** and returns `bool` indicating whether any entry was improved.
- **`getAllSubsets()` utility function** generates all non-empty subsets (bitmask enumeration). Used for both destination subsets and **forwarding set candidates** — now enumerating **all** non-empty subsets of neighbors instead of just singletons + full set (key improvement over Phase 1–4).
- **Multi-pass convergence in `solve()`**: Repeats routing table building until no entry improves (up to 50 passes), supporting undirected/cyclic graphs.
- **`TeeBuf` class** for dual-output (console + `phase5_output.txt` in append mode with run timestamps).
- **`printOptimalForwardingSets()`** method prints all forwarding sets for all nodes with costs.
- **`simulateFromForwardingToDestination()`** broadcasts to **all neighbors** of the forwarding node (wireless broadcast), with a safety limit of 500 retransmissions.
- **Phase 2 XOR recovery** now handles **multi-destination (>2)** via pairwise XOR between all destination pairs, then individual sends for remaining.

### Network Topology (Example)

```
Sources/Dests        Forwarding        Sources/Dests
  ┌───┐              ┌───┐              ┌───┐
  │ 0 │──────────────│ 2 │──────────────│ 4 │
  │ 1 │──────────────│ 3 │──────────────│ 5 │
  └───┘              └───┘              └───┘
  Sends to {4,5}                        Sends to {0,1}
```

### Input File Format (`phase5_input.txt`)

```
# Number of nodes
6
# Number of source nodes
4
# Source-Destination mapping
# Format: source_node  num_destinations  dest1 dest2 ...
0 2 4 5
1 2 4 5
4 2 0 1
5 2 0 1
# Number of edges (undirected)
8
# Edges (node1  node2  cost)
0 2 4
0 3 5.5
...
# Packets per source
10000
# Coding window size
64
```

---

## Phase 6 — Multi-Source Multi-Destination Multi-Hop (Optimal FS)

**File:** `phase6.cpp` (1134 lines) · **Input:** `phase6_input.txt`

### What's New (vs Phase 5)

- Supports **3+ hops** (multiple layers of forwarding nodes).
- Larger networks (9–12+ nodes) with multiple forwarding layers.
- **Layer-aware forwarding** in the simulation: Forwarding nodes distinguish between:
  - **Last-hop nodes** (can reach at least one destination directly via `fwdNeighbors`) → use coding window logic.
  - **Intermediate nodes** (cannot reach any destination directly) → relay packets to the next layer, behaving like a source: pop from queue, transmit to next-layer forwarding node, retransmit on failure.
- **Bottom-up routing in `solve()`** propagates routing tables from destinations backwards through **all intermediate layers** using a `while(true)` loop that repeatedly finds predecessors of the current layer, building routing tables for each new layer until no new predecessors are found.
- Uses **all non-empty subsets** of neighbors for forwarding set selection (optimal).

### Network Topology (Example — 12 nodes, 3 hops)

```
Sources (0,1,2)    Layer 1 (3,4,5)    Layer 2 (6,7,8)    Dests (9,10,11)
   ┌─┐               ┌─┐                ┌─┐               ┌──┐
   │0│───────────────│3│───────────────│6│──────────────│ 9│
   │1│───────────────│4│───────────────│7│──────────────│10│
   │2│───────────────│5│───────────────│8│──────────────│11│
   └─┘               └─┘                └─┘               └──┘
```

### Key Logic — Multi-Hop Forwarding (in simulation)

```
if (forwarding node can reach a destination directly):
    → Use CODING WINDOW logic (same as Phase 5)
else:
    → Act as RELAY: pop packet from queue, use best forwarding set
      to transmit to next-layer forwarding node. Retransmit on failure.
      If a destination receives directly during relay, mark as delivered.
```

### Input File Format (`phase6_input.txt`)

```
# Number of nodes
12
# Number of source nodes
6
# Source-Destination mapping
0 3 9 10 11
1 3 9 10 11
2 3 9 10 11
9 3 0 1 2
10 3 0 1 2
11 3 0 1 2
# Number of edges (undirected)
21
# Edges (node1  node2  cost)
0 3 4
0 4 5.5
...
# Packets per source
10000
# Coding window size
64
```

---

## Phase 6a — Multi-Hop with Heuristic (Singleton FS)

**File:** `phase6a.cpp` (1149 lines) · **Input:** `phase6a_input.txt`

### What's New (vs Phase 6)

- Introduces a **heuristic forwarding set strategy**: the routing table still stores the **optimal cost** (computed over all forwarding set subsets), but the `bestForwardingSet` used during **simulation** is restricted to the **best singleton** (single-neighbor) forwarding set.
- Adds `singletonCostTable`: `map<int, map<set<int>, double>>` — stores the cost when using the best singleton FS (for reference/comparison).
- In `buildRoutingTableForNode()`:
  - All non-empty subsets of neighbors are evaluated to compute the true optimal cost stored in `routingTable` (needed for correct downstream cost calculations).
  - Separately, the best **singleton** forwarding set (FS with `fs.size() == 1`) is tracked and stored in `bestForwardingSet` for use during simulation.
- **Purpose**: Allows comparison of simulation results between optimal (Phase 6) and heuristic-singleton (Phase 6a) forwarding strategies. The singleton heuristic is simpler and may be more practical for real implementations.

### Structural Differences from Phase 6

The only code difference is in `buildRoutingTableForNode()`:
- Phase 6 stores the overall best FS (any size) in `bestForwardingSet`.
- Phase 6a stores the best singleton FS in `bestForwardingSet`, while routing table uses the optimal cost.

---

## Phase 7 — Optimized Multi-Hop (Optimal FS + Performance)

**File:** `phase7.cpp` (1132 lines) · **Input:** `phase7_input.txt`

### What's New (vs Phase 6)

- **Performance optimization: `neighborCache`** — a pre-computed `map<int, set<int>>` neighbor adjacency map populated during `addEdge()`. The `getNeighbors()` method now returns from the cache instead of scanning all edges each time.
- **Pre-computed probability outcomes**: In the simulation, `calculateProbabilityRegions()` is called **once** before the retry loop for source-to-forwarding and forwarding-to-destination transmissions, rather than being recomputed on every retry. This is a significant performance improvement for high-packet-count simulations.
- **Increased safety limit**: `maxIter` scaled to `totalPackets * 5000` (vs `totalPackets * 500` in Phase 5/6) to support larger networks.
- **Forwarding node neighbor pre-computation**: When processing a forwarding node, `getNeighbors(selectedNode)` is called once and reused for all source iterations on that node.
- **Last-hop targeting**: `simulateFromForwardingToDestination()` is called with only the **reachable destinations** (`reachableDests` = destinations that are direct neighbors of the forwarding node), not all destinations.
- Uses **all non-empty subsets** for forwarding set selection (same optimal strategy as Phase 6).

---

## Phase 7a — Optimized Multi-Hop (Heuristic Singleton FS + Performance)

**File:** `phase7a.cpp` (1147 lines) · **Input:** `phase7a_input.txt`

### What's New (vs Phase 7)

- Combines the **performance optimizations** of Phase 7 with the **heuristic singleton forwarding set** strategy of Phase 6a.
- Contains `singletonCostTable` and uses best singleton FS in `bestForwardingSet` (same heuristic as 6a).
- Contains `neighborCache` and pre-computed probability outcomes (same optimizations as Phase 7).
- `printOptimalForwardingSets()` prints `singletonCostTable` costs instead of `routingTable` costs, showing the heuristic cost.

### Purpose

Serves as the heuristic counterpart to Phase 7 for benchmarking. Together, Phase 7 (optimal) and Phase 7a (heuristic) allow direct comparison of forwarding set strategies under the same optimized simulation engine.

---

## How to Build & Run

### Prerequisites

- C++ compiler supporting C++11 or later (e.g., `g++`, MSVC)

### Compile

```bash
# Phase 1 (interactive stdin input)
g++ -o phase1 Phase1.cpp

# Phase 2 (reads phase2_input.txt)
g++ -o phase2 phase2.cpp

# Phase 3 (interactive stdin input)
g++ -o phase3 phase3.cpp

# Phase 4 (reads phase4_input.txt)
g++ -o phase4 phase4.cpp

# Phase 5 (reads phase5_input.txt)
g++ -o phase5 phase5.cpp

# Phase 6 (reads phase6_input.txt)
g++ -o phase6 phase6.cpp

# Phase 6a (reads phase6a_input.txt)
g++ -o phase6a phase6a.cpp

# Phase 7 (reads phase7_input.txt)
g++ -o phase7 phase7.cpp

# Phase 7a (reads phase7a_input.txt)
g++ -o phase7a phase7a.cpp
```

### Run

```bash
# Phase 1 & 3 — enter values interactively
./phase1
./phase3

# Phase 2, 4, 5, 6, 6a, 7, 7a — ensure input file is in the same directory
./phase2
./phase4
./phase5
./phase6
./phase6a
./phase7
./phase7a
```

> **Note:** Phases 5–7a also save output to `phaseX_output.txt` (append mode with timestamps).

---

## Input File Format

Input files use a simple text format with `#` comment lines:

| Field | Description |
|-------|-------------|
| `# ...` | Comment line (ignored) |
| First value | Number of nodes in the network |
| Source(s) | Source node ID(s) or source-destination mapping |
| Destinations | Destination node IDs |
| Edges | `from  to  cost` — one per line |
| Packets | Number of packets to simulate per source |
| Coding window | Size of the coding window at forwarding nodes |

> **Phases 5–7a** use a source-destination mapping format:
> `source_node  num_destinations  dest1 dest2 ...`

> **Phases 5–7a** use **undirected** edges — each edge is specified once and applied in both directions with the same cost automatically.

---

## Project Structure

```
BTP/
├── Phase1.cpp            Phase 1: Core routing algorithm (stdin)
├── phase2.cpp            Phase 2: Packet simulation (file input)
├── phase2_input.txt      Sample input for Phase 2
├── phase3.cpp            Phase 3: Network coding XOR (stdin)
├── phase4.cpp            Phase 4: Multi-source NC + coding window
├── phase4_input.txt      Sample input for Phase 4
├── phase5.cpp            Phase 5: Multi-src multi-dest 2-hop (optimal FS)
├── phase5_input.txt      Sample input for Phase 5
├── phase5_output.txt     Output log for Phase 5
├── phase6.cpp            Phase 6: Multi-hop, optimal FS
├── phase6_input.txt      Sample input for Phase 6
├── phase6_output.txt     Output log for Phase 6
├── phase6a.cpp           Phase 6a: Multi-hop, heuristic singleton FS
├── phase6a_input.txt     Sample input for Phase 6a
├── phase6a_output.txt    Output log for Phase 6a
├── phase7.cpp            Phase 7: Optimized multi-hop, optimal FS
├── phase7_input.txt      Sample input for Phase 7
├── phase7_output.txt     Output log for Phase 7
├── phase7a.cpp           Phase 7a: Optimized multi-hop, heuristic singleton FS
├── phase7a_input.txt     Sample input for Phase 7a
├── phase7a_output.txt    Output log for Phase 7a
├── heuristic.txt         Reference input configurations
└── Readme.md             This file
```

---

"All rights reserved. Unauthorized use prohibited."