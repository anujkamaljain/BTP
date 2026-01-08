#include <bits/stdc++.h>
using namespace std;

// Use long double for precision
using ld = long double;
const ld INF_LD = 1e300L;   
const ld EPS = 1e-18L;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cout << "Enter number of nodes: " << flush;
    while (!(cin >> n) || n <= 0) {
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Invalid. Enter a positive integer: " << flush;
    }

    // link success probabilities p[u][v]
    vector<vector<ld>> p(n, vector<ld>(n, 0.0L));
    for (int u = 0; u < n; ++u) {
        for (int v = 0; v < n; ++v) {
            if (u == v) continue;
            char r;
            cout << "Is there a path from " << u << " -> " << v << "? (y/n): " << flush;
            cin >> r;
            if (r == 'y' || r == 'Y') {
                long double etx;
                cout << "  Enter ETX for " << u << " -> " << v << "  (ETX >= 1): " << flush;
                while (!(cin >> etx) || etx < 1.0L - (long double)EPS) {
                    cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
                    cout << "  Invalid. Enter ETX >= 1: " << flush;
                }
                ld prob = 1.0L / (ld)etx;
                prob = max<ld>(0.0L, min<ld>(1.0L, prob));
                p[u][v] = prob;
            }
        }
    }

    int src;
    cout << "\nEnter source node (0.." << n - 1 << "): " << flush;
    while (!(cin >> src) || src < 0 || src >= n) {
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Invalid. Enter source node (0.." << n - 1 << "): " << flush;
    }

    cout << "Enter number of destination nodes: " << flush;
    int m;
    while (!(cin >> m) || m <= 0 || m > 20) {
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Invalid. Enter a positive integer (<=20): " << flush;
    }

    vector<int> dests(m);
    cout << "Enter destination node IDs (space separated): " << flush;
    for (int i = 0; i < m; ++i) {
        while (!(cin >> dests[i]) || dests[i] < 0 || dests[i] >= n) {
            cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Invalid id. Enter destination node ID (0.." << n - 1 << "): " << flush;
        }
    }

    // mask of destinations per node
    vector<int> nodeDestMask(n, 0);
    for (int i = 0; i < m; ++i) nodeDestMask[dests[i]] = (1 << i);

    int FULL = (1 << m) - 1;

    // build neighbor lists
    vector<vector<int>> neigh(n);
    for (int u = 0; u < n; ++u) {
        for (int v = 0; v < n; ++v) {
            if (u != v && p[u][v] > 0.0L) neigh[u].push_back(v);
        }
    }

    // DP tables
    vector<vector<ld>> C(n, vector<ld>(1 << m, INF_LD));
    vector<vector<vector<int>>> bestFS(n, vector<vector<int>>(1 << m));
    for (int u = 0; u < n; ++u) C[u][FULL] = 0.0L;

    // process masks by descending popcount so transitions go to larger masks already computed
    vector<int> masks;
    masks.reserve(FULL + 1);
    for (int mask = 0; mask <= FULL; ++mask) masks.push_back(mask);
    sort(masks.begin(), masks.end(), [&](int a, int b) {
        int ca = __builtin_popcount((unsigned)a);
        int cb = __builtin_popcount((unsigned)b);
        return ca != cb ? ca > cb : a < b;
    });

    for (int mask : masks) {
        if (mask == FULL) continue;
        for (int u = 0; u < n; ++u) {
            int startingMask = mask;
            if ((nodeDestMask[u] & (~mask)) != 0) {
                startingMask = mask | nodeDestMask[u];
                if (startingMask == FULL) {
                    C[u][mask] = 0.0L;
                    bestFS[u][mask].clear();
                    continue;
                }
            }

            if (neigh[u].empty()) {
                C[u][mask] = INF_LD;
                bestFS[u][mask].clear();
                continue;
            }

            ld bestCost = C[u][mask];
            vector<int> bestSet;
            int d = (int)neigh[u].size();
            int maxS = 1 << d;

            for (int sMask = 1; sMask < maxS; ++sMask) {
                // prob none in S receives
                ld probNone = 1.0L;
                for (int i = 0; i < d; ++i) if (sMask & (1 << i))
                    probNone *= (1.0L - p[u][neigh[u][i]]);
                ld Psucc = 1.0L - probNone;
                if (Psucc <= EPS) continue;

                ld expectedRemain = 0.0L;
                bool bad = false;

                // iterate receiver subsets R of S
                for (int rMask = sMask; rMask; rMask = (rMask - 1) & sMask) {
                    ld probR = 1.0L;
                    int maskAfter = startingMask;
                    for (int i = 0; i < d; ++i) if (sMask & (1 << i)) {
                        int v = neigh[u][i];
                        if (rMask & (1 << i)) {
                            probR *= p[u][v];
                            maskAfter |= nodeDestMask[v];
                        } else probR *= (1.0L - p[u][v]);
                    }
                    ld cnext = C[u][maskAfter];
                    if (cnext >= INF_LD / 2) { bad = true; break; }
                    expectedRemain += probR * cnext;
                }
                if (bad) continue;

                ld candidate = (1.0L + expectedRemain) / Psucc;
                if (candidate + 1e-15L < bestCost) {
                    bestCost = candidate;
                    bestSet.clear();
                    for (int i = 0; i < d; ++i) if (sMask & (1 << i)) bestSet.push_back(neigh[u][i]);
                }
            }

            C[u][mask] = bestCost;
            bestFS[u][mask] = bestSet;
        }
    }

    // Backward sharing computation for source node
    {
        int mask = 0;
        ld bestCost = INF_LD;
        vector<int> bestSet;
        vector<int> relays;
        relays.reserve(n);

        for (int i = 0; i < n; ++i)
            if (i != src && !nodeDestMask[i] && !neigh[i].empty())
                relays.push_back(i);

        int d = (int)relays.size();
        int maxS = 1 << d;
        for (int sMask = 1; sMask < maxS; ++sMask) {
            ld probNone = 1.0L;
            for (int i = 0; i < d; ++i) if (sMask & (1 << i))
                probNone *= (1.0L - p[src][relays[i]]);
            ld Psucc = 1.0L - probNone;
            if (Psucc <= EPS) continue;

            ld expectedRemain = 0.0L;
            for (int rMask = sMask; rMask; rMask = (rMask - 1) & sMask) {
                ld probR = 1.0L;
                ld minCost = INF_LD;
                for (int i = 0; i < d; ++i) if (sMask & (1 << i)) {
                    if (rMask & (1 << i)) {
                        probR *= p[src][relays[i]];
                        minCost = min(minCost, C[relays[i]][mask]);
                    } else probR *= (1.0L - p[src][relays[i]]);
                }
                if (probR > 0) {
                    if (minCost >= INF_LD / 2) { expectedRemain = INF_LD; break; }
                    expectedRemain += probR * minCost;
                }
            }
            if (expectedRemain >= INF_LD / 2) continue;

            ld candidate = (1.0L + expectedRemain) / Psucc;
            if (candidate < bestCost) {
                bestCost = candidate;
                bestSet.clear();
                for (int i = 0; i < d; ++i) if (sMask & (1 << i)) bestSet.push_back(relays[i]);
            }
        }

        // write back so final table shows source cost
        if (bestCost < INF_LD / 2) {
            C[src][0] = bestCost;
            bestFS[src][0] = bestSet;
        }

        cout << "\nSource Backward Computation:\n";
        if (bestCost >= INF_LD / 2)
            cout << "Source " << src << " has no viable multicast path.\n";
        else {
            cout << "Minimum cost from source " << src << " to all destinations = "
                 << (double)bestCost << " expected transmissions.\n";
            cout << "Best forwarding set: { ";
            for (size_t i = 0; i < bestSet.size(); ++i)
                cout << bestSet[i] << (i + 1 < bestSet.size() ? ", " : " ");
            cout << "}\n";
        }
    }

    // Print final table
    cout << "\nRouting Table:\n\n";
    cout << left << setw(6) << "Node" << setw(16) << "Cost" << "ForwardingSet\n";
    cout << "-----------------------------------------\n";
    cout.setf(std::ios::fixed); cout << setprecision(6);
    for (int u = 0; u < n; ++u) {
        ld cost = C[u][0];
        cout << left << setw(6) << u;
        if (cost >= INF_LD / 2) cout << setw(16) << "INF";
        else cout << setw(16) << (double)cost;
        if (bestFS[u][0].empty()) {
            if (nodeDestMask[u]) cout << " (dst)";
            else cout << " -";
        } else {
            cout << " (";
            for (size_t i = 0; i < bestFS[u][0].size(); ++i) {
                cout << bestFS[u][0][i];
                if (i + 1 < bestFS[u][0].size()) cout << ",";
            }
            cout << ")";
        }
        cout << "\n";
    }

    return 0;
}