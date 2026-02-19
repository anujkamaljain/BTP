#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <climits>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <vector>

using namespace std;

const double INF = numeric_limits<double>::infinity();
const double EPSILON = 0.001;

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
    cout << "Added edge: " << from << " -> " << to << " with cost " << cost
         << endl;
  }

  void addDestination(int dest) {
    destinations.insert(dest);
    cout << "Added destination: " << dest << endl;
  }

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
  int simulatePacketTransmission(int currentNode, set<int> uncoveredDests,
                                 mt19937 &rng, int &transmissions) {
    // Base case: all destinations reached
    if (uncoveredDests.empty()) {
      return transmissions;
    }

    // Get best forwarding set for current node to uncovered destinations
    if (!bestForwardingSet[currentNode].count(uncoveredDests)) {
      // No route found - should not happen if routing table is complete
      cout << "DEBUG: No forwarding set found for node " << currentNode
           << " to {";
      for (int d : uncoveredDests)
        cout << d << " ";
      cout << "}" << endl;
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
    uniform_real_distribution<> dis(0.0, 1.0);
    double randomNum = dis(rng);

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
      return simulatePacketTransmission(currentNode, newUncovered, rng,
                                        transmissions);
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
      // No valid next node found - this happens when only destinations received
      // In this case, ret ransmit from current node to reach remaining
      // destinations
      return simulatePacketTransmission(currentNode, newUncovered, rng,
                                        transmissions);
    }

    // Continue from the best next node
    return simulatePacketTransmission(nextNode, newUncovered, rng,
                                      transmissions);
  }

  // Run simulation for multiple packets
  void runSimulation(int numPackets) {
    cout << "\n\n********** PACKET TRANSMISSION SIMULATION **********" << endl;
    cout << "Simulating " << numPackets << " packets..." << endl;

    // Initialize random number generator
    random_device rd;
    mt19937 rng(rd());

    vector<int> transmissionCounts;
    int totalTransmissions = 0;
    int minTransmissions = INT_MAX;
    int maxTransmissions = 0;

    // Simulate each packet
    for (int i = 0; i < numPackets; i++) {
      int transmissions = 0;
      int result =
          simulatePacketTransmission(source, destinations, rng, transmissions);

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
  void solve() {
    cout << "\n********** OPPORTUNISTIC ROUTING ALGORITHM **********" << endl;
    cout << "Source: " << source << endl;
    cout << "Destinations: {";
    for (int d : destinations)
      cout << d << " ";
    cout << "}" << endl;
    cout << "Number of nodes: " << numNodes << endl;

    // Step 1: Build routing tables for destination nodes
    // cout << "\n--- Step 1: Initialize destination nodes ---" << endl;
    for (int dest : destinations) {
      set<int> single;
      single.insert(dest);
      routingTable[dest][single] = 0.0;
      // cout << "Node " << dest << " is a destination, cost to itself = 0"
      //      << endl;
    }

    // Step 2: Build routing tables for neighbors of destinations
    // cout << "\n--- Step 2: Build routing tables for neighbors of destinations
    // "
    //         "---"
    //      << endl;
    set<int> predecessors = getPredecessors(destinations);

    for (int node : predecessors) {
      if (destinations.find(node) == destinations.end()) {
        buildRoutingTableForNode(node);
      }
    }

    // Step 3: Build routing table for source
    // cout << "\n--- Step 3: Build routing table for source node ---" << endl;
    buildRoutingTableForNode(source);

    cout << "\n\n********** FINAL ROUTING TABLES **********" << endl;

    for (auto &nodeEntry : routingTable) {
      int node = nodeEntry.first;
      cout << "\nNode " << node << " routing table:" << endl;

      for (auto &destEntry : nodeEntry.second) {
        const set<int> &destSet = destEntry.first;
        double cost = destEntry.second;

        cout << "  " << node << " -> {";
        for (int d : destSet)
          cout << d << " ";
        cout << "} : " << fixed << setprecision(2) << cost;

        if (bestForwardingSet[node].count(destSet)) {
          cout << " (FS={";
          for (int f : bestForwardingSet[node][destSet])
            cout << f << " ";
          cout << "})";
        }
        cout << endl;
      }
    }

    cout << "\n\n********** FINAL RESULT **********" << endl;
    cout << "Cost from source " << source << " to destination set {";
    for (int d : destinations)
      cout << d << " ";
    cout << "} = " << fixed << setprecision(2)
         << routingTable[source][destinations] << endl;
    cout << "Best Forwarding Set: {";
    for (int f : bestForwardingSet[source][destinations])
      cout << f << " ";
    cout << "}" << endl;
  }

  void displayTopology() {
    cout << "\n=== Network Topology ===" << endl;
    cout << "Node1 -> Node2 : Cost" << endl;
    for (const auto &edge : edges) {
      cout << edge.from << " -> " << edge.to << " : " << edge.cost << endl;
    }
  }
};

int main() {
  cout << "=== OPPORTUNISTIC ROUTING SIMULATION ===" << endl;
  cout << "========================================" << endl;

  int numNodes, numEdges, source, numDest;

  cout << "\nEnter number of nodes: ";
  cin >> numNodes;

  cout << "Enter source node: ";
  cin >> source;

  OpportunisticRouting router(numNodes, source);

  cout << "Enter number of destinations: ";
  cin >> numDest;

  cout << "Enter destination nodes: ";
  for (int i = 0; i < numDest; i++) {
    int dest;
    cin >> dest;
    router.addDestination(dest);
  }

  cout << "\nEnter number of edges: ";
  cin >> numEdges;

  cout << "Enter edges (from to cost):" << endl;
  for (int i = 0; i < numEdges; i++) {
    int from, to;
    double cost;
    cin >> from >> to >> cost;
    router.addEdge(from, to, cost);
  }

  router.displayTopology();
  router.solve();

  // Run packet transmission simulation
  int numPackets;
  cout << "\nEnter number of packets to simulate: ";
  cin >> numPackets;

  if (numPackets > 0) {
    router.runSimulation(numPackets);
  }

  return 0;
}
