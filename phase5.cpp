// Multisource-Multidest (4 src/dest with 2 hop) with Coding Window

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

// Generate all non-empty subsets of a set
vector<set<int>> getAllSubsets(const set<int> &s) {
  vector<int> elems(s.begin(), s.end());
  int n = elems.size();
  vector<set<int>> subsets;
  for (int mask = 1; mask < (1 << n); mask++) {
    set<int> subset;
    for (int i = 0; i < n; i++) {
      if (mask & (1 << i)) {
        subset.insert(elems[i]);
      }
    }
    subsets.push_back(subset);
  }
  return subsets;
}

class OpportunisticRouting {
private:
  int numNodes;
  vector<Edge> edges;
  map<int, map<int, double>> graph;
  map<int, map<set<int>, double>> routingTable;
  map<int, map<set<int>, set<int>>> bestForwardingSet;
  set<int> srcDestNodes;
  map<int, set<int>> sourceDestMap; // source -> its specific destinations

public:
  OpportunisticRouting(int n) : numNodes(n) {}

  void addEdge(int from, int to, double cost) {
    edges.push_back({from, to, cost});
    graph[from][to] = cost;
  }

  void addSrcDestNode(int node) { srcDestNodes.insert(node); }

  void setSourceDests(int src, const set<int>& dests) { sourceDestMap[src] = dests; }
  set<int> getSourceDests(int src) { return sourceDestMap[src]; }

  set<int> getSrcDestNodes() { return srcDestNodes; }

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
    if (destSet.size() == 1 && destSet.find(node) != destSet.end()) {
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
          double minCost = INF;

          for (int sNode : successfulNodes) {
            double costFromSNode = 0.0;

            if (uncoveredDests.size() == 1) {
              int uncovDest = *uncoveredDests.begin();
              set<int> singleDestSet;
              singleDestSet.insert(uncovDest);

              if (routingTable[sNode].count(singleDestSet)) {
                costFromSNode = routingTable[sNode][singleDestSet];
              } else if (graph[sNode].count(uncovDest)) {
                costFromSNode = graph[sNode][uncovDest];
              } else if (graph[node].count(uncovDest)) {
                costFromSNode = graph[node][uncovDest];
              } else {
                costFromSNode = INF;
              }
            } else {
              if (routingTable[sNode].count(uncoveredDests)) {
                costFromSNode = routingTable[sNode][uncoveredDests];
              } else if (routingTable[node].count(uncoveredDests)) {
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

    return x;
  }

  // Build routing table for a node for ALL subsets of destinations
  // Returns true if any routing table entry was improved
  bool buildRoutingTableForNode(int node, const set<int> &destinations) {
    bool changed = false;
    set<int> neighbors = getNeighbors(node);

    // Generate ALL non-empty subsets of the destination set
    vector<set<int>> destSubsets = getAllSubsets(destinations);

    for (const auto &destSet : destSubsets) {
      bool isDestination =
          (destSet.find(node) != destSet.end() && destSet.size() == 1);
      if (isDestination) {
        routingTable[node][destSet] = 0.0;
        continue;
      }

      double bestCost = INF;
      set<int> bestFS;

      // Enumerate ALL non-empty subsets of neighbors as forwarding set candidates
      // (singletons + full set is insufficient with undirected edges)
      vector<set<int>> forwardingSets = getAllSubsets(neighbors);

      for (const auto &fs : forwardingSets) {
        double cost = calculateCostIterative(node, destSet, fs);

        if (cost < bestCost) {
          bestCost = cost;
          bestFS = fs;
        }
      }

      // Only update if we found a better cost than existing entry
      double oldCost = routingTable[node].count(destSet) ? routingTable[node][destSet] : INF;
      if (bestCost < oldCost - EPSILON) {
        routingTable[node][destSet] = bestCost;
        bestForwardingSet[node][destSet] = bestFS;
        changed = true;
      } else if (!routingTable[node].count(destSet)) {
        // First time: store even if not improved
        routingTable[node][destSet] = bestCost;
        bestForwardingSet[node][destSet] = bestFS;
        changed = true;
      }
    }
    return changed;
  }

  // =========== SIMULATION METHODS ===========

  struct TransmissionOutcome {
    set<int> successfulNodes;
    double probability;
    double cumProbStart, cumProbEnd;
  };

  vector<TransmissionOutcome>
  calculateProbabilityRegions(int node, const set<int> &forwardingSet) {
    vector<TransmissionOutcome> outcomes;
    vector<int> fsNodes;
    vector<double> p;

    for (int neighbor : forwardingSet) {
      if (graph[node].count(neighbor)) {
        fsNodes.push_back(neighbor);
        p.push_back(1.0 / graph[node][neighbor]);
      }
    }

    int n = fsNodes.size();
    double cumProb = 0.0;

    for (int mask = 0; mask < (1 << n); mask++) {
      TransmissionOutcome outcome;
      outcome.probability = 1.0;

      for (int i = 0; i < n; i++) {
        if (mask & (1 << i)) {
          outcome.probability *= p[i];
          outcome.successfulNodes.insert(fsNodes[i]);
        } else {
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

  // Simulate packet from a source until ANY of targetDests receives it
  int simulateToAnyDestination(int sourceNode, const set<int> &targetDests, int &transmissions) {
    int currentNode = sourceNode;
    set<int> uncoveredDests = targetDests;

    while (true) {
      if (!bestForwardingSet[currentNode].count(uncoveredDests)) {
        // Try to find any subset match
        return -1;
      }

      set<int> forwardingSet = bestForwardingSet[currentNode][uncoveredDests];
      if (forwardingSet.empty())
        return -1;

      vector<TransmissionOutcome> outcomes =
          calculateProbabilityRegions(currentNode, forwardingSet);

      double randomNum = (double)randomgeneration(global_seed);

      TransmissionOutcome selectedOutcome;
      for (const auto &outcome : outcomes) {
        if (randomNum >= outcome.cumProbStart &&
            randomNum < outcome.cumProbEnd) {
          selectedOutcome = outcome;
          break;
        }
      }

      transmissions++;

      // Check if ANY target destination received
      for (int destId : targetDests) {
        if (selectedOutcome.successfulNodes.count(destId)) {
          return destId;
        }
      }

      if (selectedOutcome.successfulNodes.empty()) {
        continue;
      }

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
        continue;
      }

      currentNode = nextNode;
    }
  }

  // Simulate from forwarding node to destinations
  // Forwarding node BROADCASTS to ALL its neighbors (wireless)
  // Only destination nodes' receptions are counted
  set<int> simulateFromForwardingToDestination(int forwardingNode,
                                                const set<int> &targetDests,
                                                int &transmissions) {
    // Get ALL neighbors of the forwarding node (wireless broadcast reaches all)
    set<int> allNeighbors = getNeighbors(forwardingNode);
    if (allNeighbors.empty()) return {};

    int maxTrans = 500; // Safety limit

    while (maxTrans-- > 0) {
      vector<TransmissionOutcome> outcomes =
          calculateProbabilityRegions(forwardingNode, allNeighbors);
      double randVal = (double)randomgeneration(global_seed);

      TransmissionOutcome selected;
      for (const auto &outcome : outcomes) {
        if (randVal >= outcome.cumProbStart && randVal < outcome.cumProbEnd) {
          selected = outcome;
          break;
        }
      }

      transmissions++;

      // Only check if any DESTINATION received — ignore all other receivers
      set<int> receivedDests;
      for (int dest : targetDests) {
        if (selected.successfulNodes.count(dest)) {
          receivedDests.insert(dest);
        }
      }

      if (!receivedDests.empty()) {
        return receivedDests;
      }

      // No destination received — retry broadcast from same forwarding node
    }
    return {}; // Safety: max transmissions reached
  }

  // =========== PHASE 5: MULTI-SOURCE MULTI-DESTINATION NC SIMULATION ===========

  // Run Phase 5 Multi-Source Multi-Destination NC Simulation with Coding Window
  void runPhase5Simulation(const vector<int> &allSrcDestNodes, int packetsPerSource,
                           int codingWindowSize) {
    cout << "\n\n********** PHASE 5: MULTI-SOURCE MULTI-DESTINATION NC WITH CODING WINDOW "
            "**********"
         << endl;
    cout << "Source-Destination nodes: {";
    for (int s : allSrcDestNodes)
      cout << s << " ";
    cout << "}" << endl;
    cout << "Packets per source: " << packetsPerSource << endl;
    cout << "Coding window size: " << codingWindowSize << endl;

    // Data structures
    map<int, vector<int>> sourcePackets; // source -> [packet IDs]
    // forwardingQueues: key = pair(fwd_node, source) -> queue of (packet_id, source_id)
    map<pair<int,int>, queue<pair<int, int>>> forwardingQueues;
    // destReceived[source][dest] = {packet IDs received at dest from source}
    map<int, map<int, set<int>>> destReceived;

    // Coding window structures per forwarding node, per source
    map<int, map<int, vector<int>>> codingWindow;
    map<int, map<int, map<int, int>>> windowRecvCount;

    // Initialize source packets
    int packetId = 1;
    for (int src : allSrcDestNodes) {
      for (int i = 0; i < packetsPerSource; i++) {
        sourcePackets[src].push_back(packetId);
        packetId++;
      }
    }
    int totalPackets = packetId - 1;

    // Find forwarding (intermediate) nodes
    set<int> forwardingNodes;
    set<int> srcDestSet(allSrcDestNodes.begin(), allSrcDestNodes.end());
    for (int i = 0; i < numNodes; i++) {
      if (srcDestSet.find(i) == srcDestSet.end()) {
        bool hasEdge = false;
        for (auto &edge : edges) {
          if (edge.from == i || edge.to == i) {
            hasEdge = true;
            break;
          }
        }
        if (hasEdge) {
          forwardingNodes.insert(i);
        }
      }
    }

    cout << "Forwarding (intermediate) nodes: {";
    for (int f : forwardingNodes)
      cout << f << " ";
    cout << "}" << endl;

    // Initialize data structures
    for (int node : forwardingNodes) {
      for (int src : allSrcDestNodes) {
        forwardingQueues[{node, src}] = queue<pair<int, int>>();
        codingWindow[node][src] = vector<int>();
        for (int dest : sourceDestMap[src]) {
          windowRecvCount[node][src][dest] = 0;
        }
      }
    }

    // All selectable nodes: all nodes in the network (0,1,2,3,4,5)
    vector<int> selectableNodes;
    for (int s : allSrcDestNodes)
      selectableNodes.push_back(s);
    for (int f : forwardingNodes)
      selectableNodes.push_back(f);

    double totalTransmissions = 0.0;
    map<int, double> perSourceTransmissions;
    for (int src : allSrcDestNodes) {
      perSourceTransmissions[src] = 0.0;
    }

    // Phase 1: Queue-based transmission with Coding Window
    cout << "\nPhase 1: Queue-based transmission with Coding Window..." << endl;

    int iterCount = 0;
    int maxIter = totalPackets * 500; // Safety limit

    while (iterCount < maxIter) {
      iterCount++;

      // Check if all sources empty, all queues empty, and all coding windows empty
      bool allEmpty = true;
      for (int src : allSrcDestNodes) {
        if (!sourcePackets[src].empty()) {
          allEmpty = false;
          break;
        }
      }
      if (allEmpty) {
        for (int fwd : forwardingNodes) {
          for (int src : allSrcDestNodes) {
            if (!forwardingQueues[{fwd, src}].empty()) {
              allEmpty = false;
              break;
            }
            if (!codingWindow[fwd][src].empty()) {
              allEmpty = false;
              break;
            }
          }
          if (!allEmpty) break;
        }
      }
      if (allEmpty)
        break;

      // Random node selection from ALL selectable nodes (0,1,2,3,4,5)
      int nodeIdx = (int)(randomgeneration(global_seed) * selectableNodes.size());
      if (nodeIdx >= (int)selectableNodes.size()) nodeIdx = selectableNodes.size() - 1;
      int selectedNode = selectableNodes[nodeIdx];

      // Check if it's a source-destination node acting as source
      if (srcDestSet.count(selectedNode)) {
        if (!sourcePackets[selectedNode].empty()) {
          int pktId = sourcePackets[selectedNode].back();
          sourcePackets[selectedNode].pop_back();

          // Destinations for this source (specific mapping)
          set<int> thisDestinations = sourceDestMap[selectedNode];

          // Use best forwarding set for this source-destination pair
          set<int> neighbors;
          if (bestForwardingSet[selectedNode].count(thisDestinations)) {
            neighbors = bestForwardingSet[selectedNode][thisDestinations];
          } else {
            for (auto &edge : edges) {
              if (edge.from == selectedNode) {
                neighbors.insert(edge.to);
              }
            }
          }

          if (!neighbors.empty()) {
            int retryCount = 0;
            // Simulate transmission to neighbors
            while (retryCount < 1000) {
              retryCount++;
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
              perSourceTransmissions[selectedNode]++;

              if (!selected.successfulNodes.empty()) {
                // Check if any destination received directly
                bool destReceivedDirectly = false;
                for (int dest : thisDestinations) {
                  if (selected.successfulNodes.count(dest)) {
                    destReceived[selectedNode][dest].insert(pktId);
                    destReceivedDirectly = true;
                    break;
                  }
                }
                if (destReceivedDirectly) break;

                // No destination received, find best forwarding node
                int bestNode = -1;
                double minCost = INF;
                for (int node : selected.successfulNodes) {
                  if (forwardingNodes.count(node) &&
                      routingTable[node].count(thisDestinations)) {
                    double cost = routingTable[node][thisDestinations];
                    if (cost < minCost) {
                      minCost = cost;
                      bestNode = node;
                    }
                  }
                }
                if (bestNode != -1) {
                  forwardingQueues[{bestNode, selectedNode}].push({pktId, selectedNode});
                  break;
                }
              }
              // Retransmit if nothing useful received
            }
          }
        }
      }
      // Check if it's a forwarding node (Coding Window logic)
      else if (forwardingNodes.count(selectedNode)) {
        // Process each source's coding window independently
        for (int src : allSrcDestNodes) {
          pair<int,int> queueKey = {selectedNode, src};

          // Step 1: Fill coding window from queue
          while ((int)codingWindow[selectedNode][src].size() < codingWindowSize &&
                 !forwardingQueues[queueKey].empty()) {
            pair<int, int> windowEntry = forwardingQueues[queueKey].front();
            forwardingQueues[queueKey].pop();
            codingWindow[selectedNode][src].push_back(windowEntry.first);
          }

          // Step 2: If coding window has packets, transmit a coded packet
          if (!codingWindow[selectedNode][src].empty()) {
            int targetWindowSize = (int)codingWindow[selectedNode][src].size();

            // Destinations for this source (specific mapping)
            set<int> thisDestinations = sourceDestMap[src];

            int transmissions = 0;
            set<int> receivedDests = simulateFromForwardingToDestination(
                selectedNode, thisDestinations, transmissions);
            totalTransmissions += transmissions;
            perSourceTransmissions[src] += transmissions;



            // Increment counter for ALL destinations that received
            for (int dest : receivedDests) {
              windowRecvCount[selectedNode][src][dest]++;
            }

            // Step 3: Check if all destinations received enough coded packets
            bool allDestsReceived = true;
            for (int dest : thisDestinations) {
              if (windowRecvCount[selectedNode][src][dest] < targetWindowSize) {
                allDestsReceived = false;
                break;
              }
            }

            if (allDestsReceived && targetWindowSize > 0) {
              for (int pId : codingWindow[selectedNode][src]) {
                for (int dest : thisDestinations) {
                  destReceived[src][dest].insert(pId);
                }
              }
              codingWindow[selectedNode][src].clear();
              for (int dest : thisDestinations) {
                windowRecvCount[selectedNode][src][dest] = 0;
              }
            }

            // Step 4: Refill window
            while ((int)codingWindow[selectedNode][src].size() < codingWindowSize &&
                   !forwardingQueues[queueKey].empty()) {
              pair<int, int> newEntry = forwardingQueues[queueKey].front();
              forwardingQueues[queueKey].pop();
              codingWindow[selectedNode][src].push_back(newEntry.first);
            }
          }
        }
      }
    }

    if (iterCount >= maxIter) {
      cout << "WARNING: Reached max iterations (" << maxIter << "). Some packets may not be delivered." << endl;
    }

    double phase1Trans = totalTransmissions;
    cout << "Phase 1 transmissions: " << fixed << setprecision(2) << phase1Trans
         << endl;

    // Phase 2: XOR Recovery for each source's missing packets
    int xorCount = 0;
    int individualCount = 0;
    int edgeCount = 2;

    for (int src : allSrcDestNodes) {
      // Get destinations for this source (specific mapping)
      vector<int> destList(sourceDestMap[src].begin(), sourceDestMap[src].end());

      // Calculate packet ID range for this source
      int srcIdx = 0;
      for (int s : allSrcDestNodes) {
        if (s == src) break;
        srcIdx++;
      }
      int startPkt = srcIdx * packetsPerSource + 1;
      int endPkt = startPkt + packetsPerSource - 1;

      // For multi-destination (>2), do pairwise XOR then individual
      // First, find missing packets per destination
      map<int, set<int>> missingPerDest;
      for (int dest : destList) {
        for (int pkt = startPkt; pkt <= endPkt; pkt++) {
          if (destReceived[src][dest].count(pkt) == 0) {
            missingPerDest[dest].insert(pkt);
          }
        }
      }

      // Try pairwise XOR between destination pairs
      for (int i = 0; i < (int)destList.size(); i++) {
        for (int j = i + 1; j < (int)destList.size(); j++) {
          int d1 = destList[i];
          int d2 = destList[j];

          while (!missingPerDest[d1].empty() && !missingPerDest[d2].empty()) {
            int p1 = *missingPerDest[d1].begin();
            int p2 = *missingPerDest[d2].begin();
            totalTransmissions += edgeCount;
            perSourceTransmissions[src] += edgeCount;
            xorCount++;
            missingPerDest[d1].erase(p1);
            missingPerDest[d2].erase(p2);
            destReceived[src][d1].insert(p1);
            destReceived[src][d2].insert(p2);
          }
        }
      }

      // Individual sends for any still-missing packets
      for (int dest : destList) {
        for (int pkt : missingPerDest[dest]) {
          totalTransmissions += edgeCount;
          perSourceTransmissions[src] += edgeCount;
          individualCount++;
          destReceived[src][dest].insert(pkt);
        }
      }
    }

    // Results
    double avgTrans = totalTransmissions / packetsPerSource;

    cout << "\n********** PHASE 5 RESULTS (MULTI-SRC MULTI-DEST CODING WINDOW) **********" << endl;
    cout << "Total packets (all sources combined): " << totalPackets << endl;
    cout << "Packets per source: " << packetsPerSource << endl;
    cout << "Coding Window Size: " << codingWindowSize << endl;
    cout << "Total transmissions: " << fixed << setprecision(2)
         << totalTransmissions << endl;
    cout << "Average per packet (overall): " << fixed << setprecision(2) << avgTrans
         << endl;
  }

  // Build routing tables for ALL source-dest pairs
  void solve() {
    // Initialize routing for ALL destination nodes across all sources
    for (int src : srcDestNodes) {
      set<int> thisDestinations = sourceDestMap[src];
      for (int dest : thisDestinations) {
        set<int> single;
        single.insert(dest);
        routingTable[dest][single] = 0.0;
      }
    }

    // Multi-pass convergence: repeat until no routing table entry improves
    bool changed = true;
    int pass = 0;
    const int maxPasses = 50; // Safety limit

    while (changed && pass < maxPasses) {
      changed = false;
      pass++;

      for (int src : srcDestNodes) {
        set<int> thisDestinations = sourceDestMap[src];

        // Build routing tables for predecessors of destinations (bottom-up)
        set<int> predecessors = getPredecessors(thisDestinations);
        for (int node : predecessors) {
          if (thisDestinations.find(node) == thisDestinations.end()) {
            if (buildRoutingTableForNode(node, thisDestinations))
              changed = true;
          }
        }

        // Build routing table for the source
        if (thisDestinations.find(src) == thisDestinations.end()) {
          if (buildRoutingTableForNode(src, thisDestinations))
            changed = true;
        }
      }
    }
  }

  void printOptimalForwardingSets() {
    cout << "\n========== OPTIMAL FORWARDING SETS (All nodes) ==========" << endl;
    for (int node = 0; node < numNodes; node++) {
      if (bestForwardingSet.find(node) == bestForwardingSet.end()) continue;
      for (const auto &destEntry : bestForwardingSet[node]) {
        const set<int> &destSet = destEntry.first;
        const set<int> &fs = destEntry.second;
        cout << "Node " << node << " -> Dest{";
        for (auto it = destSet.begin(); it != destSet.end(); it++) {
          if (it != destSet.begin()) cout << ",";
          cout << *it;
        }
        cout << "} : Forwarding Set = {";
        for (auto it = fs.begin(); it != fs.end(); it++) {
          if (it != fs.begin()) cout << ",";
          cout << *it;
        }
        cout << "}";
        if (routingTable[node].count(destSet)) {
          cout << " (cost=" << fixed << setprecision(4) << routingTable[node][destSet] << ")";
        }
        cout << endl;
      }
    }
    cout << "============================================================" << endl;
  }

};

// Helper: read next non-comment, non-empty line from file
static bool getNextDataLine(ifstream &fin, string &line) {
  while (getline(fin, line)) {
    size_t start = line.find_first_not_of(" \t\r\n");
    if (start == string::npos) continue;
    line = line.substr(start);
    if (line[0] == '#') continue;
    return true;
  }
  return false;
}

int main() {
  const string inputFile = "phase5_input.txt";
  ifstream fin(inputFile);
  if (!fin.is_open()) {
    cerr << "ERROR: Could not open " << inputFile << endl;
    cerr << "Make sure phase5_input.txt is in the same directory." << endl;
    return 1;
  }

  cout << "=================================================" << endl;
  cout << "  PHASE 5: MULTI-SOURCE MULTI-DESTINATION NC" << endl;
  cout << "  All nodes act as BOTH source and destination" << endl;
  cout << "=================================================" << endl;
  cout << "Reading input from: " << inputFile << endl;

  string line;
  int numNodes, numEdges, numSrcDest;

  // Number of nodes
  getNextDataLine(fin, line);
  numNodes = stoi(line);
  cout << "Number of nodes: " << numNodes << endl;

  // Number of source-destination nodes
  getNextDataLine(fin, line);
  numSrcDest = stoi(line);
  cout << "Number of source-destination nodes: " << numSrcDest << endl;

  OpportunisticRouting router(numNodes);

  // Source-Destination mapping (each source with its specific destinations)
  vector<int> srcDestNodes;
  for (int i = 0; i < numSrcDest; i++) {
    getNextDataLine(fin, line);
    istringstream iss(line);
    int srcNode, numDests;
    iss >> srcNode >> numDests;
    set<int> dests;
    for (int j = 0; j < numDests; j++) {
      int d;
      iss >> d;
      dests.insert(d);
    }
    srcDestNodes.push_back(srcNode);
    router.addSrcDestNode(srcNode);
    router.setSourceDests(srcNode, dests);
    cout << "Source node " << srcNode << " -> Destinations: {";
    for (int d : dests) cout << d << " ";
    cout << "}" << endl;
  }

  // Number of edges
  getNextDataLine(fin, line);
  numEdges = stoi(line);
  cout << "Number of edges (undirected): " << numEdges << endl;

  // Edges — undirected, so add both directions automatically
  for (int i = 0; i < numEdges; i++) {
    getNextDataLine(fin, line);
    istringstream iss(line);
    int from, to;
    double cost;
    iss >> from >> to >> cost;
    router.addEdge(from, to, cost);  // forward
    router.addEdge(to, from, cost);  // reverse (same cost)
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

  // Build routing tables
  router.solve();
  router.printOptimalForwardingSets();
  cout.flush();

  // Run simulation
  if (packetsPerSource > 0 && codingWindowSize > 0) {
    router.runPhase5Simulation(srcDestNodes, packetsPerSource, codingWindowSize);
  }

  return 0;
}
