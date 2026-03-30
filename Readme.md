# Opportunistic Routing with Network Coding — Simulation Framework

A progressive C++ simulation of **opportunistic routing** and **network coding** techniques for wireless networks. The project is organized into 6 phases, each building upon the previous one to add more realistic and complex capabilities..

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
- [Phase 6 — Multi-Source Multi-Destination Multi-Hop](#phase-6--multi-source-multi-destination-multi-hop)
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
| **Iterative Cost Calculation** | Uses fixed-point iteration over all `2^n` reception outcomes (bitmask enumeration) to compute expected transmission cost. |
| **Coding Window** | A buffer at forwarding nodes. Packets are grouped into windows; coded (XOR'd) packets are sent. Once all destinations receive enough coded packets (`≥ window size`), they can decode all original packets. |
| **XOR Recovery** | After Phase 1 delivery, missing packets at destinations are recovered via XOR: one coded packet simultaneously satisfies one missing packet per destination. |
| **L'Ecuyer RNG** | Phases 4–6 use the L'Ecuyer combined multiplicative congruential generator for reproducible, high-quality random numbers. |

---

## Phase Summary

| Phase | File | Input | Graph | Sources | Destinations | Hops | Network Coding | Coding Window |
|-------|------|-------|-------|---------|--------------|------|----------------|---------------|
| 1 | `Phase1.cpp` | stdin | Directed | 1 | 1–2 | Any | ✗ | ✗ |
| 2 | `phase2.cpp` | `phase2_input.txt` | Directed | 1 | 1–2 | Any | ✗ | ✗ |
| 3 | `phase3.cpp` | stdin | Directed | 1 | 2 | Any | ✓ (XOR) | ✗ |
| 4 | `phase4.cpp` | `phase4_input.txt` | Directed | 1–2 | 2 | 2 | ✓ (XOR) | ✓ |
| 5 | `phase5.cpp` | `phase5_input.txt` | **Undirected** | 2–4 | 2+ each | 2 | ✓ (XOR) | ✓ |
| 6 | `phase6.cpp` | `phase6_input.txt` | **Undirected** | 3–6 | 3 each | **3** | ✓ (XOR) | ✓ |

---

## Phase 1 — Core Routing Algorithm

**File:** `Phase1.cpp` · **Input:** stdin (interactive prompts)

### What It Does

Implements the foundational opportunistic routing algorithm. Given a directed network graph, a source, and a set of destinations, it computes:

1. **Routing tables** — expected transmission cost from every relevant node to every destination subset.
2. **Best forwarding sets** — the optimal set of neighbors to forward to for minimum cost.

### Key Logic

- **Cost computation** (`calculateCostIterative`): For each candidate forwarding set, enumerates all `2^n` possible reception outcomes (success/failure per neighbor) using bitmask iteration, computes the probability-weighted expected transmissions, and converges via fixed-point iteration (tolerance `ε = 0.001`, max 1000 iterations).
- **Bottom-up table building** (`solve`): First initializes destination nodes (cost = 0), then computes costs for predecessors of destinations, and finally for the source.
- Prints the **final routing tables** and **best forwarding set** for the source.

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

### Output

- Network topology display
- Per-node routing tables with forwarding sets
- Final cost and best forwarding set from source to destinations

---

## Phase 2 — Packet Transmission Simulation

**File:** `phase2.cpp` · **Input:** `phase2_input.txt`

### What's New (vs Phase 1)

- Reads input from a **text file** (with `#` comment support) instead of stdin.
- Adds a **Monte Carlo packet simulation** that sends thousands of packets through the network to **empirically validate** the theoretical routing costs.

### Key Logic

- **Probability regions** (`calculateProbabilityRegions`): Maps the `[0, 1)` random number space to each of the `2^n` reception outcomes using cumulative probability ranges.
- **Packet simulation** (`simulatePacketTransmission`): For each packet, generates a random number, determines which neighbors received, removes covered destinations, and recursively continues from the best receiving node until all destinations are reached.
- Reports **avg / min / max transmissions** across all packets and compares to the theoretical cost.

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
- Theoretical vs empirical cost comparison

---

## Phase 3 — Network Coding (XOR)

**File:** `phase3.cpp` · **Input:** stdin (interactive)

### What's New (vs Phase 2)

- Introduces **network coding simulation** using XOR-based recovery.
- Adds a `simulateToAnyDestination` method — packet is forwarded until **any one** destination receives it (not all).

### Key Logic — Network Coding Simulation

The simulation operates in **two phases per batch**:

| Sub-Phase | Description |
|-----------|-------------|
| **Phase 1** | Each packet is routed until **any one** destination receives it. Track which destination got which packet. |
| **Phase 2 (Recovery)** | Find missing packets at each destination. Use **XOR recovery**: if dest₁ has packet A (missing at dest₂) and dest₂ has packet B (missing at dest₁), then `XOR(A, B)` satisfies both in a single coded transmission. Remaining gaps are filled individually. |

> Requires exactly **2 destinations**.

### Input Format (stdin)

Same as Phase 1, plus prompts for:
- Number of batches
- Packets per batch

### Output

- Total packets processed
- XOR transmissions used
- Individual transmissions used
- Average transmissions per packet **with** network coding

---

## Phase 4 — Multi-Source NC with Coding Window

**File:** `phase4.cpp` · **Input:** `phase4_input.txt`

### What's New (vs Phase 3)

- **Multiple sources** (1–2) feeding packets into the network.
- **Coding window** mechanism at forwarding (intermediate) nodes.
- Replaces `std::mt19937` with the **L'Ecuyer combined multiplicative congruential RNG** for higher-quality simulation randomness.
- File-based input with separate source and destination specifications.

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

1. **Source → Forwarding:** Sources pop packets and transmit probabilistically. Forwarding nodes queue received packets.
2. **Coding Window:** At each forwarding node, packets fill a window (up to `W` size). The node transmits **coded** packets to destinations. Each coded packet is a random linear combination.
3. **Window Flush:** Once **every destination** has received `≥ W` coded packets, all `W` original packets are considered decoded/delivered. The window is cleared and refilled.
4. **XOR Recovery (Phase 2):** After all windows are flushed, any remaining missing packets are handled via pairwise XOR or individual retransmission.

### Node Selection

A random node is selected each iteration from `sources ∪ forwarding_nodes`. If it's a source, it sends a new packet. If it's a forwarding node, it processes its coding window.

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

### Output

- Phase 1 total transmissions
- Phase 2 XOR and individual recovery counts
- Average transmissions per packet

---

## Phase 5 — Multi-Source Multi-Destination (2-Hop)

**File:** `phase5.cpp` · **Input:** `phase5_input.txt`

### What's New (vs Phase 4)

- Nodes act as **both source and destination** simultaneously (e.g., node 0 sends to {4,5} while node 4 sends to {0,1}).
- **Undirected graph** — edges are added in both directions with the same cost.
- Each source has its **own destination set** (source-destination mapping).
- Coding windows and forwarding queues are maintained **per (forwarding_node, source)** pair.
- Routing tables are computed for **all non-empty subsets** of each source's destinations using `getAllSubsets()`.
- Supports 4+ source-destination nodes.

### Network Topology (Example)

```
Sources/Dests        Forwarding        Sources/Dests
  ┌───┐              ┌───┐              ┌───┐
  │ 0 │──────────────│ 2 │──────────────│ 4 │
  │ 1 │──────────────│ 3 │──────────────│ 5 │
  └───┘              └───┘              └───┘
  Sends to {4,5}                        Sends to {0,1}
```

### Key Logic

- **Solve:** For each source, builds routing tables bottom-up from its destinations through predecessors to the source itself.
- **Simulation:** Randomly selects a node. If it's a source-dest node with packets, transmits towards its destinations. If it's a forwarding node, processes each source's coding window independently.
- **Wireless broadcast:** Forwarding nodes broadcast to **all** neighbors; only destination receptions are tracked.

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

### Output

- Optimal forwarding sets for each source
- Phase 1 transmissions
- XOR recovery and individual sends per source
- Total and average per-packet transmissions

---

## Phase 6 — Multi-Source Multi-Destination Multi-Hop

**File:** `phase6.cpp` · **Input:** `phase6_input.txt`

### What's New (vs Phase 5)

- Supports **3+ hops** (multiple layers of forwarding nodes).
- Larger networks (9–12+ nodes) with 3 forwarding layers.
- **Layer-aware forwarding:** Forwarding nodes distinguish between:
  - **Last-hop nodes** (can reach destinations directly) → use coding window.
  - **Intermediate nodes** (relay packets to the next layer) → forward packets to the next closer forwarding node, similar to source transmission.
- **Bottom-up routing:** `solve()` propagates routing tables from destinations backwards through **all intermediate layers** until the source is reached.

### Network Topology (Example — 12 nodes, 3 hops)

```
Sources (0,1,2)    Layer 1 (3,4,5)    Layer 2 (6,7,8)    Dests (9,10,11)
   ┌─┐               ┌─┐                ┌─┐               ┌──┐
   │0│───────────────│3│───────────────│6│──────────────│ 9│
   │1│───────────────│4│───────────────│7│──────────────│10│
   │2│───────────────│5│───────────────│8│──────────────│11│
   └─┘               └─┘                └─┘               └──┘
```

### Key Logic — Multi-Hop Forwarding

```
if (forwarding node can reach a destination directly):
    → Use CODING WINDOW logic (same as Phase 5)
else:
    → Act as RELAY: pop packet from queue, transmit to next-layer
      forwarding node or directly to destination if possible
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

### Output

- Optimal forwarding sets for all source nodes
- Phase 1 transmissions (queue-based with coding window)
- Phase 2 XOR and individual recovery
- Total and average per-packet transmissions

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
```

### Run

```bash
# Phase 1 & 3 — enter values interactively
./phase1
./phase3

# Phase 2, 4, 5, 6 — ensure input file is in the same directory
./phase2
./phase4
./phase5
./phase6
```

---

## Input File Format

Input files use a simple text format with `#` comment lines:

| Field | Description |
|-------|-------------|
| `# ...` | Comment line (ignored) |
| First value | Number of nodes in the network |
| Source(s) | Source node ID(s) |
| Destinations | Destination node IDs |
| Edges | `from  to  cost` — one per line |
| Packets | Number of packets to simulate per source |
| Coding window | Size of the coding window at forwarding nodes |

> **Phases 5 & 6** use a source-destination mapping format:  
> `source_node  num_destinations  dest1 dest2 ...`

> **Phases 5 & 6** use **undirected** edges — each edge is specified once and applies in both directions with the same cost.

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
├── phase5.cpp            Phase 5: Multi-src multi-dest 2-hop
├── phase5_input.txt      Sample input for Phase 5
├── phase6.cpp            Phase 6: Multi-src multi-dest multi-hop
├── phase6_input.txt      Sample input for Phase 6
├── heuristic.txt         Reference input configurations
└── Readme.md             This file
```

---

"All rights reserved. Unauthorized use prohibited."