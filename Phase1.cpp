#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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

  return 0;
}