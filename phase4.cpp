//Network Coding

#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

const double INF = numeric_limits<double>::infinity();
const double EPSILON = 0.001;

// Required constants/macros for L'Ecuyer RNG
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS_RNG 1.2e-7
#define RNMX (1.0-EPS_RNG)

float randomgeneration(long &idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (idum <= 0)
    {
        if (-(idum) < 1)
            idum = 1;
        else
            idum = -(idum);

        idum2 = idum;

        for (j = NTAB + 7; j >= 0; j--)
        {
            k = idum / IQ1;
            idum = IA1 * (idum - k * IQ1) - k * IR1;
            if (idum < 0) idum += IM1;
            if (j < NTAB) iv[j] = idum;
        }
        iy = iv[0];
    }

    k = idum / IQ1;
    idum = IA1 * (idum - k * IQ1) - k * IR1;
    if (idum < 0) idum += IM1;

    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0) idum2 += IM2;

    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;

    if (iy < 1) iy += IMM1;

    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

// Global seed - negative value triggers initialization on first call
long global_seed = -(long)time(NULL);

struct Edge {
  int from, to;
  double cost;
};

class OpportunisticRouting {
private:
  int numNodes;
  int source;
  set<int> destinations;
  vector<Edge> edges;
  map<int, map<int, double>> graph;
  map<int, map<set<int>, double>> routingTable;
  map<int, map<set<int>, set<int>>> bestForwardingSet;

public:
  OpportunisticRouting(int n, int src) : numNodes(n), source(src) {}

  void addEdge(int from, int to, double cost) {
    edges.push_back({from, to, cost});
    graph[from][to] = cost;
  }

  void addDestination(int dest) { destinations.insert(dest); }

  set<int> getNeighbors(int node) {
    set<int> neighbors;
    for (const auto &edge : edges) {
      if (edge.from == node) {
        neighbors.insert(edge.to);
      }
    }
    return neighbors;
  }

  set<int> getPredecessors(const set<int> &destSet) {
    set<int> predecessors;
    for (const auto &edge : edges) {
      if (destSet.find(edge.to) != destSet.end()) {
        predecessors.insert(edge.from);
      }
    }
    return predecessors;
  }

  double calculateCostIterative(int node, const set<int> &destSet,
                                const set<int> &forwardingSet) {
    // cout << "\n  Calculating cost from node " << node << " to dest set {";
    // for (int d : destSet)
    //   cout << d << " ";
    // cout << "} using FS {";
    // for (int f : forwardingSet)
    //   cout << f << " ";
    // cout << "}" << endl;

    if (destSet.size() == 1 && destSet.find(node) != destSet.end()) {
      // cout << "  Node " << node << " is the only destination, cost = 0"
      //      << endl;
      return 0.0;
    }

    vector<double> costs;
    vector<int> fsNodes;

    for (int neighbor : forwardingSet) {
      if (graph[node].count(neighbor)) {
        costs.push_back(graph[node][neighbor]);
        fsNodes.push_back(neighbor);
      }
    }

    if (costs.empty()) {
      return INF;
    }

    vector<double> p;
    for (double c : costs) {
      p.push_back(1.0 / c);
    }

    // cout << "  Raw probabilities (1/cost): ";
    // for (int i = 0; i < fsNodes.size(); i++) {
    //   cout << "p" << (i + 1) << "(" << fsNodes[i] << ")=" << fixed
    //        << setprecision(3) << p[i] << " ";
    // }
    // cout << endl;

    double x = 5.0;
    double prevX;
    int iterations = 0;
    const int maxIterations = 1000;

    do {
      prevX = x;
      double newX = 0.0;

      int n = fsNodes.size();
      for (int mask = 0; mask < (1 << n); mask++) {
        double probability = 1.0;
        set<int> successfulNodes;

        for (int i = 0; i < n; i++) {
          if (mask & (1 << i)) {
            probability *= p[i];
            successfulNodes.insert(fsNodes[i]);
          } else {
            probability *= (1.0 - p[i]);
          }
        }

        double remainingCost = 0.0;

        set<int> uncoveredDests = destSet;
        for (int sNode : successfulNodes) {
          uncoveredDests.erase(sNode);
        }

        if (uncoveredDests.empty()) {
          remainingCost = 0.0;
        } else if (successfulNodes.empty()) {
          remainingCost = prevX;
        } else {
          // Some nodes received, find min cost from successful nodes to
          // uncovered destinations
          double minCost = INF;

          for (int sNode : successfulNodes) {
            // Cost from this successful forwarding node to reach uncovered
            // destinations
            double costFromSNode = 0.0;

            if (uncoveredDests.size() == 1) {
              // Single uncovered destination
              int uncovDest = *uncoveredDests.begin();
              set<int> singleDestSet;
              singleDestSet.insert(uncovDest);

              // Check routing table of successful node
              if (routingTable[sNode].count(singleDestSet)) {
                costFromSNode = routingTable[sNode][singleDestSet];
              } else if (graph[sNode].count(uncovDest)) {
                // Direct link from forwarding node
                costFromSNode = graph[sNode][uncovDest];
              } else if (graph[node].count(uncovDest)) {
                // No route from forwarding node, use direct link from original
                // node
                costFromSNode = graph[node][uncovDest];
              } else {
                costFromSNode = INF;
              }
            } else {
              // Multiple uncovered destinations
              if (routingTable[sNode].count(uncoveredDests)) {
                costFromSNode = routingTable[sNode][uncoveredDests];
              } else if (routingTable[node].count(uncoveredDests)) {
                // No route from forwarding node, use routing table from
                // original node
                costFromSNode = routingTable[node][uncoveredDests];
              } else {
                costFromSNode = INF;
              }
            }

            minCost = min(minCost, costFromSNode);
          }

          remainingCost = minCost;
        }

        newX += probability * (1.0 + remainingCost);
      }

      x = newX;
      iterations++;

    } while (fabs(x - prevX) > EPSILON && iterations < maxIterations);

    // cout << "  Converged after " << iterations
    //      << " iterations, cost = " << fixed << setprecision(2) << x << endl;
    return x;
  }

  void buildRoutingTableForNode(int node) {
    // cout << "\n=== Building routing table for node " << node << " ==="
    //      << endl;

    set<int> neighbors = getNeighbors(node);
    vector<set<int>> destSubsets;

    for (int dest : destinations) {
      set<int> single;
      single.insert(dest);
      destSubsets.push_back(single);
    }

    if (destinations.size() > 1) {
      destSubsets.push_back(destinations);
    }

    for (const auto &destSet : destSubsets) {
      // cout << "\nDestination set: {";
      // for (int d : destSet)
      //   cout << d << " ";
      // cout << "}" << endl;

      bool isDestination =
          (destSet.find(node) != destSet.end() && destSet.size() == 1);
      if (isDestination) {
        routingTable[node][destSet] = 0.0;
        // cout << "Node " << node << " is the destination, cost = 0" << endl;
        continue;
      }

      double bestCost = INF;
      set<int> bestFS;

      vector<set<int>> forwardingSets;

      for (int neighbor : neighbors) {
        set<int> fs;
        fs.insert(neighbor);
        forwardingSets.push_back(fs);
      }

      if (neighbors.size() > 1) {
        forwardingSets.push_back(neighbors);
      }

      for (const auto &fs : forwardingSets) {
        double cost = calculateCostIterative(node, destSet, fs);

        if (cost < bestCost) {
          bestCost = cost;
          bestFS = fs;
        }
      }

      routingTable[node][destSet] = bestCost;
      bestForwardingSet[node][destSet] = bestFS;

      // cout << "Best FS for node " << node << " to {";
      // for (int d : destSet)
      //   cout << d << " ";
      // cout << "} is {";
      // for (int f : bestFS)
      //   cout << f << " ";
      // cout << "} with cost " << fixed << setprecision(2) << bestCost << endl;
    }
  }

  // =========== SIMULATION METHODS ===========

  struct TransmissionOutcome {
    set<int> successfulNodes;
    double probability;
    double cumProbStart, cumProbEnd;
  };

  // Calculate probability regions for random number mapping
  vector<TransmissionOutcome>
  calculateProbabilityRegions(int node, const set<int> &forwardingSet) {
    vector<TransmissionOutcome> outcomes;
    vector<int> fsNodes;
    vector<double> p;

    // Get forwarding set nodes and their probabilities
    for (int neighbor : forwardingSet) {
      if (graph[node].count(neighbor)) {
        fsNodes.push_back(neighbor);
        p.push_back(1.0 / graph[node][neighbor]);
      }
    }

    int n = fsNodes.size();
    double cumProb = 0.0;

    // Generate all 2^n outcomes
    for (int mask = 0; mask < (1 << n); mask++) {
      TransmissionOutcome outcome;
      outcome.probability = 1.0;

      for (int i = 0; i < n; i++) {
        if (mask & (1 << i)) {
          // Success
          outcome.probability *= p[i];
          outcome.successfulNodes.insert(fsNodes[i]);
        } else {
          // Failure
          outcome.probability *= (1.0 - p[i]);
        }
      }

      outcome.cumProbStart = cumProb;
      cumProb += outcome.probability;
      outcome.cumProbEnd = cumProb;
      outcomes.push_back(outcome);
    }

    return outcomes;
  }

  // Simulate transmission of one packet from currentNode to uncoveredDests
  int simulatePacketTransmission(int currentNode, set<int> uncoveredDests, int &transmissions) {
    // Base case: all destinations reached
    if (uncoveredDests.empty()) {
      return transmissions;
    }

    // Get best forwarding set for current node to uncovered destinations
    if (!bestForwardingSet[currentNode].count(uncoveredDests)) {
      // No route found - should not happen if routing table is complete
      return -1;
    }

    set<int> forwardingSet = bestForwardingSet[currentNode][uncoveredDests];

    if (forwardingSet.empty()) {
      // No forwarding set - should not happen
      return -1;
    }

    // Calculate probability regions
    vector<TransmissionOutcome> outcomes =
        calculateProbabilityRegions(currentNode, forwardingSet);

    // Generate random number [0, 1]
    double randomNum = (double)randomgeneration(global_seed);

    // Find which outcome occurred
    TransmissionOutcome selectedOutcome;
    for (const auto &outcome : outcomes) {
      if (randomNum >= outcome.cumProbStart && randomNum < outcome.cumProbEnd) {
        selectedOutcome = outcome;
        break;
      }
    }

    // Increment transmission counter
    transmissions++;

    // Update uncovered destinations (remove any that were reached and are
    // destinations)
    set<int> newUncovered = uncoveredDests;
    for (int node : selectedOutcome.successfulNodes) {
      newUncovered.erase(
          node); // Remove if it's a destination (cumulative tracking)
    }

    // Check if we're done
    if (newUncovered.empty()) {
      return transmissions;
    }

    // If no nodes received, retransmit from same node
    if (selectedOutcome.successfulNodes.empty()) {
      return simulatePacketTransmission(currentNode, newUncovered, transmissions);
    }

    // Some nodes received - find the one with minimum cost to remaining
    // destinations
    int nextNode = -1;
    double minCost = INF;

    for (int sNode : selectedOutcome.successfulNodes) {
      double cost = INF;

      // Check if this node has a route to remaining destinations
      if (routingTable[sNode].count(newUncovered)) {
        cost = routingTable[sNode][newUncovered];
      } else if (newUncovered.size() == 1) {
        // Try single destination
        set<int> singleDest;
        singleDest.insert(*newUncovered.begin());
        if (routingTable[sNode].count(singleDest)) {
          cost = routingTable[sNode][singleDest];
        }
      }

      if (cost < minCost) {
        minCost = cost;
        nextNode = sNode;
      }
    }

    if (nextNode == -1) {
      // No valid next node found - retransmit from current node
      return simulatePacketTransmission(currentNode, newUncovered, transmissions);
    }

    // Continue from the best next node
    return simulatePacketTransmission(nextNode, newUncovered, transmissions);
  }

  // Simulate packet until ANY destination receives it
  // Returns: destination ID that received, -1 on error
  int simulateToAnyDestination(int &transmissions) {
    int currentNode = source;
    set<int> uncoveredDests = destinations;

    while (true) {
      // Get forwarding set
      if (!bestForwardingSet[currentNode].count(uncoveredDests)) {
        return -1; // No route
      }

      set<int> forwardingSet = bestForwardingSet[currentNode][uncoveredDests];
      if (forwardingSet.empty())
        return -1;

      // Calculate probability regions
      vector<TransmissionOutcome> outcomes =
          calculateProbabilityRegions(currentNode, forwardingSet);

      // Generate random number
      double randomNum = (double)randomgeneration(global_seed);

      // Find outcome
      TransmissionOutcome selectedOutcome;
      for (const auto &outcome : outcomes) {
        if (randomNum >= outcome.cumProbStart &&
            randomNum < outcome.cumProbEnd) {
          selectedOutcome = outcome;
          break;
        }
      }

      transmissions++;

      // Check if ANY destination received - STOP here!
      for (int destId : destinations) {
        if (selectedOutcome.successfulNodes.count(destId)) {
          return destId; // Return which destination received
        }
      }

      // No destination received - check if forwarding node received
      if (selectedOutcome.successfulNodes.empty()) {
        // Retransmit from same node
        continue;
      }

      // Find best forwarding node to continue
      int nextNode = -1;
      double minCost = INF;
      for (int sNode : selectedOutcome.successfulNodes) {
        if (routingTable[sNode].count(uncoveredDests)) {
          double cost = routingTable[sNode][uncoveredDests];
          if (cost < minCost) {
            minCost = cost;
            nextNode = sNode;
          }
        }
      }

      if (nextNode == -1) {
        // Retransmit from current node
        continue;
      }

      currentNode = nextNode;
    }
  }

  // Run simulation for multiple packets
  void runSimulation(int numPackets) {
    cout << "\n\n********** PACKET TRANSMISSION SIMULATION **********" << endl;
    cout << "Simulating " << numPackets << " packets..." << endl;

    vector<int> transmissionCounts;
    int totalTransmissions = 0;
    int minTransmissions = INT_MAX;
    int maxTransmissions = 0;

    // Simulate each packet
    for (int i = 0; i < numPackets; i++) {
      int transmissions = 0;
      int result =
          simulatePacketTransmission(source, destinations, transmissions);

      if (result == -1) {
        cout << "Error: Simulation failed for packet " << (i + 1) << endl;
        return;
      }

      transmissionCounts.push_back(transmissions);
      totalTransmissions += transmissions;
      minTransmissions = min(minTransmissions, transmissions);
      maxTransmissions = max(maxTransmissions, transmissions);
    }

    // Calculate average
    double avgTransmissions = (double)totalTransmissions / numPackets;

    // Display results
    cout << "\n********** SIMULATION RESULTS **********" << endl;
    cout << "Number of packets simulated: " << numPackets << endl;
    cout << "Average transmissions per packet: " << fixed << setprecision(2)
         << avgTransmissions << endl;
    cout << "Min transmissions: " << minTransmissions << endl;
    cout << "Max transmissions: " << maxTransmissions << endl;
    cout << "\nTheoretical cost (from routing table): " << fixed
         << setprecision(2) << routingTable[source][destinations] << endl;
    cout << "Difference: " << fixed << setprecision(2)
         << abs(avgTransmissions - routingTable[source][destinations]) << endl;
  }

  // =========== SIMPLIFIED NETWORK CODING SIMULATION ===========

  // Run simplified network coding simulation
  void runNetworkCodingSimulation(int numBatches, int packetsPerBatch) {
    cout << "\n\n********** NETWORK CODING SIMULATION **********" << endl;
    cout << "Batches: " << numBatches
         << ", Packets per batch: " << packetsPerBatch << endl;

    if (destinations.size() != 2) {
      cout << "Network coding requires exactly 2 destinations." << endl;
      return;
    }

    vector<int> destList(destinations.begin(), destinations.end());
    int dest1 = destList[0];
    int dest2 = destList[1];

    // Get theoretical cost to both destinations (for display only)
    double costToBoth = routingTable[source][destinations];

    double totalTransmissionsNC = 0.0;
    int xorCount = 0;
    int individualCount = 0;

    // Process each batch
    for (int batch = 0; batch < numBatches; batch++) {
      map<int, set<int>>
          receivedPackets; // receivedPackets[dest_id] = {packet_ids}

      // Phase 1: Simulate each packet until ANY destination receives
      for (int pkt = 1; pkt <= packetsPerBatch; pkt++) {
        int transmissions = 0;
        int receivedBy = simulateToAnyDestination(transmissions);
        totalTransmissionsNC += transmissions;

        if (receivedBy != -1) {
          receivedPackets[receivedBy].insert(pkt);
        }
      }

      // Phase 2: Find missing packets at each destination
      set<int> dest1Missing, dest2Missing;
      for (int pkt = 1; pkt <= packetsPerBatch; pkt++) {
        if (receivedPackets[dest1].count(pkt) == 0)
          dest1Missing.insert(pkt);
        if (receivedPackets[dest2].count(pkt) == 0)
          dest2Missing.insert(pkt);
      }

      // Phase 2: Use XOR or individual transmissions (using theoretical costs)
      while (!dest1Missing.empty() || !dest2Missing.empty()) {
        bool foundXOR = false;

        // Look for XOR opportunity
        if (!dest1Missing.empty() && !dest2Missing.empty()) {
          for (int p1 : receivedPackets[dest1]) {
            if (dest2Missing.count(p1)) {
              for (int p2 : receivedPackets[dest2]) {
                if (dest1Missing.count(p2)) {
                  // XOR(p1, p2) helps both!
                  // Use cost to reach one node (XOR packet goes once, helps
                  // both)
                  // Edge count from source to destination (e.g., 0→1→3 = 2
                  // edges)
                  int edgeCount =
                      2; // Number of hops from source to destinations
                  totalTransmissionsNC += edgeCount;
                  xorCount++;

                  // Both destinations now have both packets
                  dest1Missing.erase(p2);
                  dest2Missing.erase(p1);
                  receivedPackets[dest1].insert(p2);
                  receivedPackets[dest2].insert(p1);
                  foundXOR = true;
                  break;
                }
              }
              if (foundXOR)
                break;
            }
          }
        }

        if (!foundXOR) {
          // No XOR opportunity - send missing packets individually
          if (!dest1Missing.empty()) {
            int pkt = *dest1Missing.begin();
            int edgeCount = 2; // Number of hops
            totalTransmissionsNC += edgeCount;
            individualCount++;
            dest1Missing.erase(pkt);
            receivedPackets[dest1].insert(pkt);
          } else if (!dest2Missing.empty()) {
            int pkt = *dest2Missing.begin();
            int edgeCount = 2; // Number of hops
            totalTransmissionsNC += edgeCount;
            individualCount++;
            dest2Missing.erase(pkt);
            receivedPackets[dest2].insert(pkt);
          }
        }
      }
    }

    int totalPackets = numBatches * packetsPerBatch;
    double avgNC = totalTransmissionsNC / totalPackets;

    cout << "\n********** RESULTS **********" << endl;
    cout << "Total packets: " << totalPackets << endl;
    cout << "Theoretical cost (to both dests): " << fixed << setprecision(2)
         << costToBoth << endl;
    cout << "XOR transmissions used: " << xorCount << endl;
    cout << "Individual transmissions used: " << individualCount << endl;
    cout << "Average transmissions per packet (WITH NC): " << fixed
         << setprecision(2) << avgNC << endl;
  }

  // =========== PHASE 4: MULTI-SOURCE NETWORK CODING SIMULATION ===========

  // Simulate from source to any forwarding node, return node that received
  int simulateToForwardingNode(int sourceNode, int &transmissions) {
    set<int> forwardingNodes;
    for (auto &edge : edges) {
      if (edge.from == sourceNode) {
        forwardingNodes.insert(edge.to);
      }
    }

    if (forwardingNodes.empty())
      return -1;

    // Use ETX-based transmission
    while (true) {
      vector<TransmissionOutcome> outcomes =
          calculateProbabilityRegions(sourceNode, forwardingNodes);
      double randVal = (double)randomgeneration(global_seed);

      TransmissionOutcome selected;
      for (const auto &outcome : outcomes) {
        if (randVal >= outcome.cumProbStart && randVal < outcome.cumProbEnd) {
          selected = outcome;
          break;
        }
      }

      transmissions++;

      if (!selected.successfulNodes.empty()) {
        // Find best forwarding node (min cost to destinations)
        int bestNode = -1;
        double minCost = INF;
        for (int node : selected.successfulNodes) {
          if (routingTable[node].count(destinations)) {
            double cost = routingTable[node][destinations];
            if (cost < minCost) {
              minCost = cost;
              bestNode = node;
            }
          }
        }
        if (bestNode != -1)
          return bestNode;
      }
      // Retransmit if nothing received
    }
  }

  // Simulate from forwarding node to destinations
  // Returns set of ALL destinations that received the coded packet
  set<int> simulateFromForwardingToDestination(int forwardingNode, int &transmissions) {
    int currentNode = forwardingNode;

    while (true) {
      if (!bestForwardingSet[currentNode].count(destinations)) {
        return {};  // empty set
      }

      set<int> fwdSet = bestForwardingSet[currentNode][destinations];
      if (fwdSet.empty())
        return {};

      vector<TransmissionOutcome> outcomes =
          calculateProbabilityRegions(currentNode, fwdSet);
      double randVal = (double)randomgeneration(global_seed);

      TransmissionOutcome selected;
      for (const auto &outcome : outcomes) {
        if (randVal >= outcome.cumProbStart && randVal < outcome.cumProbEnd) {
          selected = outcome;
          break;
        }
      }

      transmissions++;

      // Check ALL destinations that received this transmission
      set<int> receivedDests;
      for (int dest : destinations) {
        if (selected.successfulNodes.count(dest)) {
          receivedDests.insert(dest);
        }
      }

      // If any destination received, return ALL that received
      if (!receivedDests.empty()) {
        return receivedDests;
      }

      // No destination received - find best node to continue
      if (!selected.successfulNodes.empty()) {
        int nextNode = -1;
        double minCost = INF;
        for (int node : selected.successfulNodes) {
          if (routingTable[node].count(destinations)) {
            double cost = routingTable[node][destinations];
            if (cost < minCost) {
              minCost = cost;
              nextNode = node;
            }
          }
        }
        if (nextNode != -1) {
          currentNode = nextNode;
        }
      }
      // If no progress, continue from same node (retransmit)
    }
  }

  // Run Phase 4 Multi-Source NC Simulation with Coding Window
  void runPhase4Simulation(const vector<int> &sources, int packetsPerSource,
                           int codingWindowSize) {
    cout << "\n\n********** PHASE 4: MULTI-SOURCE NC WITH CODING WINDOW "
            "**********"
         << endl;
    cout << "Sources: {";
    for (int s : sources)
      cout << s << " ";
    cout << "}, Destinations: {";
    for (int d : destinations)
      cout << d << " ";
    cout << "}" << endl;
    cout << "Packets per source: " << packetsPerSource << endl;

    // Data structures
    map<int, vector<int>> sourcePackets; // source -> [packet IDs]
    map<int, queue<pair<int, int>>>
        forwardingQueues; // node -> queue of (packet_id, source_id)
    map<int, set<int>> destReceived; // dest -> {packet IDs received}

    // Coding window structures
    map<int, vector<int>>
        codingWindow; // forwarding node -> packet IDs in current window
    map<int, map<int, int>>
        windowRecvCount; // forwarding node -> dest -> coded packets received

    // Initialize source packets
    int packetId = 1;
    for (int src : sources) {
      for (int i = 0; i < packetsPerSource; i++) {
        sourcePackets[src].push_back(packetId++);
      }
    }
    int totalPackets = packetId - 1;

    // Find forwarding nodes (nodes that are not sources or destinations)
    set<int> forwardingNodes;
    for (auto &edge : edges) {
      if (find(sources.begin(), sources.end(), edge.to) == sources.end() &&
          destinations.find(edge.to) == destinations.end()) {
        forwardingNodes.insert(edge.to);
      }
    }

    // Initialize forwarding queues and coding windows
    for (int node : forwardingNodes) {
      forwardingQueues[node] = queue<pair<int, int>>();
      codingWindow[node] = vector<int>();
      for (int dest : destinations) {
        windowRecvCount[node][dest] = 0;
      }
    }

    // All selectable nodes (sources + forwarding)
    vector<int> selectableNodes;
    for (int s : sources)
      selectableNodes.push_back(s);
    for (int f : forwardingNodes)
      selectableNodes.push_back(f);

    double totalTransmissions = 0.0;

    // Phase 1: Queue-based transmission with Coding Window
    cout << "\nPhase 1: Queue-based transmission with Coding Window..."
         << endl;

    while (true) {
      // Check if all sources empty, all queues empty, and all coding windows
      // empty
      bool allEmpty = true;
      for (int src : sources) {
        if (!sourcePackets[src].empty())
          allEmpty = false;
      }
      for (int fwd : forwardingNodes) {
        if (!forwardingQueues[fwd].empty())
          allEmpty = false;
        if (!codingWindow[fwd].empty())
          allEmpty = false;
      }
      if (allEmpty)
        break;

      // Random node selection from sources + forwarding nodes
      int nodeIdx = (int)(randomgeneration(global_seed) * selectableNodes.size());
      if (nodeIdx >= (int)selectableNodes.size()) nodeIdx = selectableNodes.size() - 1;
      int selectedNode = selectableNodes[nodeIdx];

      // Check if it's a source node
      if (find(sources.begin(), sources.end(), selectedNode) != sources.end()) {
        if (!sourcePackets[selectedNode].empty()) {
          int pktId = sourcePackets[selectedNode].back();
          sourcePackets[selectedNode].pop_back();

          // Get all neighbors of source
          set<int> neighbors;
          for (auto &edge : edges) {
            if (edge.from == selectedNode) {
              neighbors.insert(edge.to);
            }
          }

          if (!neighbors.empty()) {
            // Simulate transmission to neighbors
            while (true) {
              vector<TransmissionOutcome> outcomes =
                  calculateProbabilityRegions(selectedNode, neighbors);
              double randVal = (double)randomgeneration(global_seed);

              TransmissionOutcome selected;
              for (const auto &outcome : outcomes) {
                if (randVal >= outcome.cumProbStart &&
                    randVal < outcome.cumProbEnd) {
                  selected = outcome;
                  break;
                }
              }

              totalTransmissions++;

              if (!selected.successfulNodes.empty()) {
                // Check if any destination received directly
                for (int dest : destinations) {
                  if (selected.successfulNodes.count(dest)) {
                    destReceived[dest].insert(pktId);
                    goto next_iteration; // Packet delivered to destination
                  }
                }

                // No destination received, find best forwarding node
                int bestNode = -1;
                double minCost = INF;
                for (int node : selected.successfulNodes) {
                  if (forwardingNodes.count(node) &&
                      routingTable[node].count(destinations)) {
                    double cost = routingTable[node][destinations];
                    if (cost < minCost) {
                      minCost = cost;
                      bestNode = node;
                    }
                  }
                }
                if (bestNode != -1) {
                  forwardingQueues[bestNode].push({pktId, selectedNode});
                  goto next_iteration;
                }
              }
              // Retransmit if nothing useful received
            }
          }
        }
      }
      // Check if it's a forwarding node (Coding Window logic)
      else if (forwardingNodes.count(selectedNode)) {
        // Step 1: Fill coding window from queue (up to codingWindowSize)
        while ((int)codingWindow[selectedNode].size() < codingWindowSize &&
               !forwardingQueues[selectedNode].empty()) {
          pair<int, int> windowEntry = forwardingQueues[selectedNode].front();
          forwardingQueues[selectedNode].pop();
          codingWindow[selectedNode].push_back(windowEntry.first);
        }

        // Step 2: If coding window has packets, transmit a coded packet
        if (!codingWindow[selectedNode].empty()) {
          // Capture window size BEFORE transmission as the decoding target
          int targetWindowSize = (int)codingWindow[selectedNode].size();

          int transmissions = 0;
          set<int> receivedDests = simulateFromForwardingToDestination(
              selectedNode, transmissions);
          totalTransmissions += transmissions;

          // Increment counter for ALL destinations that received
          for (int dest : receivedDests) {
            windowRecvCount[selectedNode][dest]++;
          }

          // Step 3: Check if all destinations received enough coded packets
          // Use the target size captured before transmission (not after refill)
          bool allDestsReceived = true;
          for (int dest : destinations) {
            if (windowRecvCount[selectedNode][dest] < targetWindowSize) {
              allDestsReceived = false;
              break;
            }
          }

          if (allDestsReceived && targetWindowSize > 0) {
            // All destinations can decode all packets in the window
            // Mark all window packets as delivered to all destinations
            for (int pId : codingWindow[selectedNode]) {
              for (int dest : destinations) {
                destReceived[dest].insert(pId);
              }
            }
            // Empty the coding window and reset counts
            codingWindow[selectedNode].clear();
            for (int dest : destinations) {
              windowRecvCount[selectedNode][dest] = 0;
            }
          }

          // Step 4: Refill window with new arrivals AFTER flush check
          while ((int)codingWindow[selectedNode].size() < codingWindowSize &&
                 !forwardingQueues[selectedNode].empty()) {
            pair<int, int> newEntry = forwardingQueues[selectedNode].front();
            forwardingQueues[selectedNode].pop();
            codingWindow[selectedNode].push_back(newEntry.first);
          }
        }
      }

    next_iteration:;
    }

    double phase1Trans = totalTransmissions;
    cout << "Phase 1 transmissions: " << fixed << setprecision(2) << phase1Trans
         << endl;

    // Phase 2: XOR Recovery
    cout << "\nPhase 2: XOR Recovery..." << endl;
    int xorCount = 0;
    int individualCount = 0;
    int edgeCount = 2; // Edges from forwarding to destinations

    vector<int> destList(destinations.begin(), destinations.end());
    if (destList.size() == 2) {
      int dest1 = destList[0];
      int dest2 = destList[1];

      // Find missing packets at each destination
      set<int> dest1Missing, dest2Missing;
      for (int pkt = 1; pkt <= totalPackets; pkt++) {
        if (destReceived[dest1].count(pkt) == 0)
          dest1Missing.insert(pkt);
        if (destReceived[dest2].count(pkt) == 0)
          dest2Missing.insert(pkt);
      }

      // XOR: each XOR delivers one unique packet to each destination
      while (!dest1Missing.empty() && !dest2Missing.empty()) {
        int p1 = *dest1Missing.begin();
        int p2 = *dest2Missing.begin();
        totalTransmissions += edgeCount;
        xorCount++;
        dest1Missing.erase(p1);
        dest2Missing.erase(p2);
        destReceived[dest1].insert(p1);
        destReceived[dest2].insert(p2);
      }

      // Individual sends for remaining
      for (int pkt : dest1Missing) {
        totalTransmissions += edgeCount;
        individualCount++;
      }
      for (int pkt : dest2Missing) {
        totalTransmissions += edgeCount;
        individualCount++;
      }
    } else if (destList.size() == 1) {
      // Single destination - individual sends for missing
      int dest = destList[0];
      for (int pkt = 1; pkt <= totalPackets; pkt++) {
        if (destReceived[dest].count(pkt) == 0) {
          totalTransmissions += edgeCount;
          individualCount++;
        }
      }
    }

    // Results
    double avgTrans = totalTransmissions / totalPackets;

    cout << "\n********** PHASE 4 RESULTS (CODING WINDOW) **********" << endl;
    cout << "Total packets: " << totalPackets << endl;
    cout << "Coding Window Size: " << codingWindowSize << endl;
    cout << "Phase 1 transmissions: " << fixed << setprecision(2) << phase1Trans
         << endl;
    cout << "Phase 2 XOR transmissions: " << xorCount << " (cost "
         << xorCount * edgeCount << ")" << endl;
    cout << "Phase 2 Individual transmissions: " << individualCount << " (cost "
         << individualCount * edgeCount << ")" << endl;
    cout << "Total transmissions: " << fixed << setprecision(2)
         << totalTransmissions << endl;
    cout << "Average per packet: " << fixed << setprecision(2) << avgTrans
         << endl;
  }

  void solve() {

    // Step 1: Build routing tables for destination nodes
    for (int dest : destinations) {
      set<int> single;
      single.insert(dest);
      routingTable[dest][single] = 0.0;
    }

    // Step 2: Build routing tables for neighbors of destinations
    set<int> predecessors = getPredecessors(destinations);
    for (int node : predecessors) {
      if (destinations.find(node) == destinations.end()) {
        buildRoutingTableForNode(node);
      }
    }

    // Step 3: Build routing table for source
    buildRoutingTableForNode(source);
  }

  void displayTopology() {
    cout << "\n=== Network Topology ===" << endl;
    cout << "Node1 -> Node2 : Cost" << endl;
    for (const auto &edge : edges) {
      cout << edge.from << " -> " << edge.to << " : " << edge.cost << endl;
    }
  }
};

// Helper: read next non-comment, non-empty line from file
static bool getNextDataLine(ifstream &fin, string &line) {
  while (getline(fin, line)) {
    // Trim leading whitespace
    size_t start = line.find_first_not_of(" \t\r\n");
    if (start == string::npos) continue; // blank line
    line = line.substr(start);
    if (line[0] == '#') continue; // comment line
    return true;
  }
  return false;
}

int main() {
  const string inputFile = "phase4_input.txt";
  ifstream fin(inputFile);
  if (!fin.is_open()) {
    cerr << "ERROR: Could not open " << inputFile << endl;
    cerr << "Make sure phase4_input.txt is in the same directory." << endl;
    return 1;
  }

  cout << "=== PHASE 4: MULTI-SOURCE NETWORK CODING ==" << endl;
  cout << "Reading input from: " << inputFile << endl;
  cout << "=============================================" << endl;

  string line;
  int numNodes, numEdges, numSources, numDest;

  // Number of nodes
  getNextDataLine(fin, line);
  numNodes = stoi(line);
  cout << "Number of nodes: " << numNodes << endl;

  // Number of sources
  getNextDataLine(fin, line);
  numSources = stoi(line);
  cout << "Number of sources: " << numSources << endl;

  // Source nodes
  vector<int> sources;
  getNextDataLine(fin, line);
  {
    istringstream iss(line);
    for (int i = 0; i < numSources; i++) {
      int src;
      iss >> src;
      sources.push_back(src);
      cout << "Source node: " << src << endl;
    }
  }

  // Use first source for routing table calculation
  OpportunisticRouting router(numNodes, sources[0]);

  // Number of destinations
  getNextDataLine(fin, line);
  numDest = stoi(line);
  cout << "Number of destinations: " << numDest << endl;

  // Destination nodes
  getNextDataLine(fin, line);
  {
    istringstream iss(line);
    for (int i = 0; i < numDest; i++) {
      int dest;
      iss >> dest;
      router.addDestination(dest);
      cout << "Destination node: " << dest << endl;
    }
  }

  // Number of edges
  getNextDataLine(fin, line);
  numEdges = stoi(line);
  cout << "Number of edges: " << numEdges << endl;

  // Edges
  for (int i = 0; i < numEdges; i++) {
    getNextDataLine(fin, line);
    istringstream iss(line);
    int from, to;
    double cost;
    iss >> from >> to >> cost;
    router.addEdge(from, to, cost);
  }

  // Packets per source
  int packetsPerSource = 0;
  if (getNextDataLine(fin, line)) {
    packetsPerSource = stoi(line);
  }

  // Coding window size
  int codingWindowSize = 0;
  if (getNextDataLine(fin, line)) {
    codingWindowSize = stoi(line);
  }

  fin.close();

  router.solve();

  if (packetsPerSource > 0 && codingWindowSize > 0) {
    cout << "\nPackets per source: " << packetsPerSource << endl;
    cout << "Coding window size: " << codingWindowSize << endl;
    router.runPhase4Simulation(sources, packetsPerSource, codingWindowSize);
  }

  return 0;
}