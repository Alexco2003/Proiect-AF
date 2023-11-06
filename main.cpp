#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <set>
#include <string>
#include <windows.h>

using namespace std;

void clearScreen()
{
    system("CLS");
}

void displayMenu()
{
    cout << "---User Menu---" << endl << endl;
    cout << "1-886. Possible Bipartition" << endl;
    cout << "2-210. Course Schedule II" << endl;
    cout << "3-1192. Critical Connections in a Network" << endl;
    cout << "4-802. Find Eventual Safe States" << endl;
    cout << "5-934. Shortest Bridge" << endl;
    cout << "6-990. Satisfiability of Equality Equations" << endl;
    cout << "7-Padure" << endl;
    cout << "8-Patrol2" << endl;
    cout << "9-F. Graph Without Long Directed Paths" << endl;
    cout << "10-F. Minimum Maximum Distance" << endl;
    cout << "11-Exit" << endl;
}

class Graph
{
private:

    int n; /// The number of nodes of the graph.
    vector<vector<int>> adjList; /// This vector stores the adjacency list representation of the graph.

    /// Private member function to check if the graph is bipartite using BFS (Breadth-First Search)
    bool isBipartite(vector<int>& colors);

    /// Private member function to perform Depth-First Search (DFS) and find topological order starting from a given node.
    bool dfsTopological(int node, vector<int>& visited, stack<int>& s);

    ///Private DFS function to find critical connections in the graph.
    void dfsCriticalConnections(int node, int parent, int& time, vector<int>& discoveryTime, vector<int>& lowTime, vector<vector<int>>& result);

    /// Private DFS function to check if a node is eventually safe.
    bool dfsSafeNodes(int node, unordered_map<int,bool>& safe);

    /// Private helper function to check if the given coordinates are invalid (out of grid boundaries).
    bool invalid(int n, int r, int c);

    /// Private DFS function to mark all cells of the first island as visited.
    void dfsBridge(int n, int r, int c, set<pair<int, int>>& visit, vector<vector<int>>& grid, vector<vector<int>>& direct);

    /// Private BFS function to find the shortest bridge between the islands.
    void bfsBridge(int n, set<pair<int, int>>& visit, vector<vector<int>>& grid, vector<vector<int>>& direct, int& res);

    /// Private DFS to check if there's a path between 'node' and 'target'
    bool dfsEquations(int node, int target, vector<bool>& visited);

    /// Private helper function to check if a block is inside the forest grid.
    bool insideForest(int N, int M, pair<int, int> block);

    /// Private BFS function to find the lowest cost path in the forest
    int bfsForest(vector<vector<int>>& forest, int N, int M, pair<int, int> start, pair<int, int> end, vector<vector<int>>& direct);

    /// Private BFS function starting from the given node and update the distances in the 'dist' matrix.
    void bfsPatrol2(int node, vector<vector<int>>& dist);

    /// Private BFS traversal function to calculate distances from given starting node to all other nodes.
    vector<int> bfsMinimumMaximum(int node);



public:

    /// Constructor to initialize the graph with given edges and directed flag.
    Graph(int n, vector<vector<int>>& edges, bool directed);
    /// Constructor without parameters.
    Graph();

    /// Public member function to check if it's possible to split the graph into two groups.
    /// https://leetcode.com/problems/possible-bipartition/
    bool possibleBipartition();

    /// Public member function to find the topological order of nodes in the graph.
    /// \return A vector representing the topological order of nodes, or an empty vector if the graph contains a cycle.
    /// https://leetcode.com/problems/course-schedule-ii/
    vector<int> findOrder();

    /// Public function to find critical connections in the graph.
    /// \return Vector containing critical connections (bridges) in the graph.
    /// https://leetcode.com/problems/critical-connections-in-a-network
    vector<vector<int>> criticalConnections();

    /// Public function to find all nodes that are eventually safe in the graph.
    /// \return A vector containing all nodes that are eventually safe.
    /// https://leetcode.com/problems/find-eventual-safe-states/
    vector<int> eventualSafeNodes();

    /// Public function to find the shortest bridge connecting the two islands in the grid.
    /// \param grid Input grid representing the matrix.
    /// \return The minimum number of flips needed to connect the two islands.
    /// https://leetcode.com/problems/shortest-bridge/
    int shortestBridge(vector<vector<int>>& grid);

    /// Public function to check if given equations are possible.
    /// https://leetcode.com/problems/satisfiability-of-equality-equations/
    bool equationsPossible(vector<string>& equations);

    /// Public function to find the lowest cost path between the prince and the castle.
    /// \return The lowest cost path
    /// https://www.infoarena.ro/problema/padure
    int padure(int N, int M, int pl, int pc, int cl, int cc, vector<vector<int>>& forest);

    /// Public function to determine the orientation of edges in a graph to avoid paths of length two or greater.
    /// \param orderList A vector containing the order of edges as they were given in the input.
    /// https://codeforces.com/contest/1144/problem/F
    void orientation(vector<int>& orderList);

    /// Public function to find the minimum time required to escape considering all patrol constraints.
    /// \param dist: Reference to the 2D vector representing distances between nodes at different times.
    /// \return The minimum time required to escape or -1 if not possible.
    /// https://www.infoarena.ro/problema/patrol2
    int patrol2(vector<vector<int>>& dist);

    /// Public function to calculate the minimum maximum distance among marked nodes in the graph.
    /// \param k The number of marked nodes.
    /// \param marked A boolean vector indicating whether a node is marked or not.
    /// \return The minimum maximum distance among the marked nodes.
    /// https://codeforces.com/contest/1881/problem/F
    int minimumMaximum(int k, vector<bool>& marked);

};


/// Constructor to initialize the graph with given edges and directed flag.
Graph::Graph(int n, vector<vector<int>>& edges, bool directed)
{
    this->n = n; /// Initialize the number of nodes.
    this->adjList.resize(this->n + 1); /// Resize the adjacency list to accommodate 'n' nodes.

    /// Iterate through the given edges and populate the adjacency list.
    for (const auto& edge : edges)
    {
        this->adjList[edge[0]].push_back(edge[1]); /// Add edge from node 'edge[0]' to node 'edge[1]'.

        /// If the graph is undirected, add reverse edge from 'edge[1]' to 'edge[0]'.
        if (directed == false)
        {
            this->adjList[edge[1]].push_back(edge[0]);
        }
    }
}

/// Constructor without parameters.
Graph::Graph()
{
    this->n = 0; /// Initialize the number of nodes to 0 for the default constructor.
    this->adjList = vector<vector<int>>(); /// Initialize the adjacency list as an empty vector of vectors.
}

/// Private member function to check if the graph is bipartite using BFS (Breadth-First Search)
/// \param colors Colors array
bool Graph::isBipartite(vector<int>& colors)
{
    queue<int> q; /// Queue to perform BFS traversal of the graph.
    for (int i = 1; i <= this->n; i++)
    {
        if (colors[i] != -1) /// If the color is already assigned, skip this node.
            continue;
        q.push(i); /// Push the current node into the queue to start BFS from this node.
        colors[i] = 0; /// Assign color 0 to the current node.

        while (!q.empty()) /// Perform BFS traversal until the queue is empty.
        {
            int curr = q.front(); /// Get the front node from the queue.
            q.pop(); /// Remove the front node from the queue.

            for (int neighbor : this->adjList[curr]) /// Iterate through neighbors of the current node.
            {
                if (colors[neighbor] == -1) /// If neighbor's color is not assigned, assign the opposite color to it.
                {
                    colors[neighbor] = 1 - colors[curr]; /// Assign opposite color to the neighbor.
                    q.push(neighbor); /// Push the neighbor into the queue for further exploration.
                }
                else
                {
                    if (colors[neighbor] == colors[curr]) /// If neighbor has the same color, the graph is not bipartite.
                        return false;
                }
            }
        }
    }
    return true; /// If all nodes can be colored without conflicts, the graph is bipartite.
}


/// Private member function to perform Depth-First Search (DFS) and find topological order starting from a given node.
/// \param node The current node being explored.
/// \param visited A vector indicating whether each node has been visited (0: not visited, 1: visiting, 2: visited).
/// \param s A stack to store the topological order of nodes.
/// \return True if the graph is acyclic and a valid topological order is found, false otherwise.
bool Graph::dfsTopological(int node, vector<int>& visited, stack<int>& s)
{
    visited[node] = 1; /// Mark the current node as being visited (visiting).

    /// Explore neighbors of the current node.
    for (int neighbor : this->adjList[node])
    {
        if (visited[neighbor] == 1)
        {
            return false; /// If the neighbor is currently being visited, a cycle is detected, and the graph is not acyclic.
        }
        if (visited[neighbor] == 0 && !dfsTopological(neighbor, visited, s))
        {
            return false; /// If the neighbor has not been visited yet and a cycle is detected in its subtree, return false.
        }
    }

    visited[node] = 2; /// Mark the current node as visited.

    s.push(node); /// Push the current node onto the stack to record its topological order.

    return true; /// The current node and its subtree have been explored without detecting a cycle.
}


///Private DFS function to find critical connections in the graph.
/// \param node Current node being explored.
/// \param parent Parent node of the current node in the DFS traversal.
/// \param time Current time during DFS traversal.
/// \param discoveryTime Vector to store discovery time of each node.
/// \param lowTime Vector to store the lowest discovery time reachable from the node.
/// \param result Vector to store critical connections.
void Graph::dfsCriticalConnections(int node, int parent, int& time, vector<int>& discoveryTime, vector<int>& lowTime, vector<vector<int>>& result)
{
    discoveryTime[node] = time;  /// Set the discovery time of the current node.
    lowTime[node] = time;  /// Set the lowest discovery time of the current node.
    time++;  /// Increment the time for the next node.

    for (int neighbor : this->adjList[node])
    {
        if (neighbor == parent)
        {
            continue; /// Skip the parent node in the DFS traversal.
        }

        if (discoveryTime[neighbor] == -1)
        {
            dfsCriticalConnections(neighbor, node, time, discoveryTime, lowTime, result);

            /// Update lowTime of the current node based on its child.
            lowTime[node] = min(lowTime[node], lowTime[neighbor]);

            /// If the edge is a bridge, add it to the result.
            if (lowTime[neighbor] > discoveryTime[node])
            {
                result.push_back({node, neighbor});
            }
        }
        else
        {
            /// Update lowTime of the current node based on the visited neighbor.
            lowTime[node] = min(lowTime[node], discoveryTime[neighbor]);
        }
    }
}


/// Private DFS function to check if a node is eventually safe.
/// \param node The current node being checked for safety.
/// \param safe A map to store the safe status of nodes.
/// \return True if the node is eventually safe, false otherwise.
bool Graph::dfsSafeNodes(int node, unordered_map<int,bool>& safe)
{
    /// If the safe status of the current node is already calculated, return it.
    if (safe.find(node) != safe.end())
        return safe[node];

    /// Mark the current node as unsafe initially.
    safe[node] = false;

    /// Explore neighbors of the current node.
    for (int neighbor : this->adjList[node])
        /// If any neighbor is not safe, the current node is not safe either.
        if (!dfsSafeNodes(neighbor, safe))
            return false;

    /// If all neighbors are safe, mark the current node as safe.
    safe[node] = true;
    return true;
}


/// Private helper function to check if the given coordinates are invalid (out of grid boundaries).
/// \param n Size of the grid.
/// \param r Row coordinate.
/// \param c Column coordinate.
/// \return true if the coordinates are invalid, false otherwise.
bool Graph::invalid(int n, int r, int c)
{
    return r < 0 || c < 0 || r >= n || c >= n;
}

/// Private DFS function to mark all cells of the first island as visited.
/// \param n Size of the grid.
/// \param r Row coordinate.
/// \param c Column coordinate.
/// \param visit Set to store visited cells.
/// \param grid Input grid representing the matrix.
/// \param direct Vector containing possible directions for exploration.
void Graph::dfsBridge(int n, int r, int c, set<pair<int, int>>& visit, vector<vector<int>>& grid, vector<vector<int>>& direct)
{
    /// Base cases: if coordinates are out of grid or cell is not land (1) or already visited
    if (invalid(n, r, c) || grid[r][c] != 1 || visit.find({r, c}) != visit.end())
    {
        return;
    }

    /// Mark the current cell as visited
    visit.insert({r, c});

    /// Explore neighbors in all directions
    for (const auto& dir : direct)
    {
        int newR = r + dir[0];
        int newC = c + dir[1];
        dfsBridge(n, newR, newC, visit, grid, direct);
    }
}

/// Private BFS function to find the shortest bridge between the islands.
/// \param n Size of the grid.
/// \param visit Set containing visited cells.
/// \param grid Input grid representing the matrix.
/// \param direct Vector containing possible directions for exploration.
/// \param res Reference variable to store the minimum number of flips needed to connect the islands.
void Graph::bfsBridge(int n, set<pair<int, int>>& visit, vector<vector<int>>& grid, vector<vector<int>>& direct, int& res)
{
    queue<pair<int, int>> q;

    /// Add all initially visited cells to the queue
    for (const auto& point : visit)
    {
        q.push(point);
    }

    while (!q.empty())
    {
        int qSize = q.size();

        /// Process nodes at the current level of BFS
        for (int i = 0; i < qSize; i++)
        {
            int r = q.front().first;
            int c = q.front().second;

            q.pop();

            /// Explore neighbors in all directions
            for (const auto& dir : direct)
            {
                int newR = r + dir[0];
                int newC = c + dir[1];

                /// If neighbor is valid and not visited, mark it as visited and enqueue for the next level of BFS
                if (!invalid(n, newR, newC) && visit.find({newR, newC}) == visit.end())
                {
                    if (grid[newR][newC] == 1)
                    {
                        return; /// The other island is found, exit
                    }
                    q.push({newR, newC});
                    visit.insert({newR, newC});
                }
            }
        }

        res++; /// Increment the number of flips needed
    }
}

/// Private DFS to check if there's a path between 'node' and 'target'
/// \param node Current node in the DFS traversal.
/// \param target The target node to find a path to.
/// \param visited A boolean array to mark visited nodes.
/// \return True if a path is found from 'node' to 'target', false otherwise.
bool Graph::dfsEquations(int node, int target, vector<bool>& visited)
{
    /// If the node equals the target, a path is found
    if (node == target) return true;

    /// Mark the current node as visited
    visited[node] = true;

    /// Explore neighbors of the current node
    for (int neighbor : this->adjList[node])
    {
        /// If the neighbor is not visited, explore it
        if (!visited[neighbor])
        {
            if (dfsEquations(neighbor, target, visited))
            {
                return true; /// If a path is found, return true
            }
        }
    }
    return false; /// If no path is found, return false
}


/// Private helper function to check if a block is inside the forest grid.
/// \param N Number of rows in the forest grid.
/// \param M Number of columns in the forest grid.
/// \param block Coordinates of the block to be checked.
/// \return true if the block is inside the grid, false otherwise.
bool Graph::insideForest(int N, int M, pair<int, int> block)
{
    return (block.first >= 1 && block.second >= 1 && block.first <= N && block.second <= M);
}

/// Private BFS function to find the lowest cost path in the forest
/// \param forest Forest grid represented as a matrix.
/// \param N Number of rows in the forest grid.
/// \param M Number of columns in the forest grid.
/// \param start Starting position coordinates.
/// \param end Destination position coordinates.
/// \param direct Possible directions for movement (up, down, left, right).
/// \return The lowest cost path from the start to the end position.
int Graph::bfsForest(vector<vector<int>>& forest, int N, int M, pair<int, int> start, pair<int, int> end, vector<vector<int>>& direct)
{
    vector<vector<int>> costs(N + 1, vector<int>(M + 1, INT_MAX)); /// Matrix of costs between start and another block
    deque<pair<int, int>> q; /// Queue for BFS traversal

    q.push_back(start); /// Push in the queue the starting position
    costs[start.first][start.second] = 0; /// Cost of the starting position is 0

    while (!q.empty())
    {
        int r = q.front().first; /// Current Row
        int c = q.front().second; /// Current Col
        q.pop_front();

        /// Explore neighbors
        for (const auto& dir : direct)
        {
            int newR = r + dir[0]; /// New Row after the move
            int newC = c + dir[1]; /// New Column after the move

            if (insideForest(N, M, {newR, newC})) /// If the new position is inside the forest grid
            {
                /// Calculate the cost of moving from the current position to the new position.
                /// If the tree type is the same, the cost is 0 (no additional cost). Otherwise, the cost is 1.
                int cost = (forest[r][c] == forest[newR][newC]) ? 0 : 1;

                /// Calculate the total cost from the start position to the current position plus the cost of the move.
                int totalCost = costs[r][c] + cost;

                /// If the newly calculated total cost is smaller than the previous cost to reach this position,
                /// Update the cost matrix with the lower cost, and add the new position to the queue for further exploration.
                if (totalCost < costs[newR][newC])
                {
                    /// Update the cost matrix with the lower total cost to reach this position.
                    costs[newR][newC] = totalCost;

                    /// If the cost of the move is 0 (zero-cost edge), prioritize it by adding the new position to the front of the queue.
                    if (cost == 0)
                    {
                        q.push_front({newR, newC});
                    }
                    /// If the cost of the move is 1 (non-zero-cost edge), add the new position to the back of the queue.
                    else
                    {
                        q.push_back({newR, newC});
                    }
                }
            }
        }
    }

    return costs[end.first][end.second]; /// Return the lowest cost path to the destination block
}

/// Private BFS function starting from the given node and update the distances in the 'dist' matrix.
/// \param node The starting node for BFS traversal.
/// \param dist The 2D vector representing distances between nodes at different times.
void Graph::bfsPatrol2(int node, vector<vector<int>>& dist)
{
    /// Create a queue to store nodes and their corresponding time
    queue<pair<int, int>> q;
    /// Push the starting node and time 0 into the queue
    q.push(make_pair(node, 0));
    /// Initialize the distance to the starting node at time 0 as 0
    dist[node][0] = 0;

    /// BFS until the queue is empty
    while (!q.empty())
    {
        /// Get the current node and time from the front of the queue
        int currNode = q.front().first;
        int currTime = q.front().second;
        q.pop();

        /// Iterate through neighbors of the current node
        for (auto neighbor : this->adjList[currNode])
        {
            /// Calculate the new time considering the cyclic nature (MCM)
            int newTime = (currTime + 1) % 420;

            /// If the neighbor node is available at the new time
            if (dist[neighbor][newTime] != -1)
            {
                /// If the new time is shorter, update the distance and enqueue the neighbor
                if (dist[neighbor][newTime] > dist[currNode][currTime] + 1)
                {
                    dist[neighbor][newTime] = dist[currNode][currTime] + 1;
                    q.push(make_pair(neighbor, newTime));
                }
            }
        }
    }
}

/// Private BFS traversal function to calculate distances from given starting node to all other nodes.
/// \param node The starting node for BFS traversal.
/// \return A vector containing distances from the starting node to all other nodes.
vector<int> Graph::bfsMinimumMaximum(int node)
{
    queue<int> q;
    q.push(node); /// Push the starting node into the queue.

    vector<int> dist(this->n + 1, -1); /// Initialize a vector to store distances from the starting node to all the others, initially set to -1 .
    dist[node] = 0; /// The distance from the starting node to itself is 0.

    while (!q.empty()) /// Perform BFS traversal to calculate distances from the starting node to all other nodes.
    {
        int curr = q.front(); /// Get the front element from the queue.
        q.pop(); /// Remove the front element from the queue.

        for (int neighbor : this->adjList[curr]) /// Iterate through neighbors of the current node.
        {
            if (dist[neighbor] == -1) /// If the neighbor node is not visited yet.
            {
                dist[neighbor] = dist[curr] + 1; /// Update the distance to the neighbor node.
                q.push(neighbor); /// Push the neighbor node into the queue for further exploration.
            }
        }
    }
    return dist; /// Return the vector of distances from the starting node to all other nodes.
}




/// Public member function to check if it's possible to split the graph into two groups
/// https://leetcode.com/problems/possible-bipartition/
bool Graph::possibleBipartition()
{
    vector<int> colors(this->n + 1, -1); /// Initialize colors array with -1 (indicating color not assigned).
    return isBipartite(colors); /// Call the private function to check if the graph is bipartite.
}

/// Public member function to find the topological order of nodes in the graph.
/// \return A vector representing the topological order of nodes, or an empty vector if the graph contains a cycle.
/// https://leetcode.com/problems/course-schedule-ii/
vector<int> Graph::findOrder()
{
    vector<int> result; /// Vector to store the resulting topological order.
    vector<int> visited(this->n, 0); /// Vector to track visited status of nodes during DFS.
    stack<int> s; /// Stack to store nodes in topological order.

    /// Iterate through each node in the graph and perform DFS to find topological order.
    for (int i = 0; i < this->n; i++)
    {
        if (visited[i] == 0 && !dfsTopological(i, visited, s))
        {
            return {}; /// If a cycle is detected during DFS, return an empty vector to indicate impossibility of topological order.
        }
    }

    /// Pop nodes from the stack and construct the topological order vector.
    while (!s.empty())
    {
        result.push_back(s.top()); /// Add the top node from the stack to the result vector.
        s.pop(); /// Remove the top node from the stack.
    }

    return result; /// Return the topological order of nodes in the graph.
}

/// Public function to find critical connections in the graph.
/// \return Vector containing critical connections (bridges) in the graph.
/// https://leetcode.com/problems/critical-connections-in-a-network
vector<vector<int>> Graph::criticalConnections()
{
    vector<int> discoveryTime(n, -1);  /// Store the discovery time of each node.
    vector<int> lowTime(n, -1);  /// Store the lowest discovery time reachable from the node.
    vector<vector<int>> result;  /// Store critical connections.

    int time = 0;  /// Variable to keep track of the current time during DFS.

    /// Perform DFS from each node to find critical connections.
    for (int i = 0; i < this->n; i++)
    {
        if (discoveryTime[i] == -1)
        {
            dfsCriticalConnections(i, -1, time, discoveryTime, lowTime, result);
        }
    }

    return result;  /// Return the vector containing critical connections.
}

/// Public function to find all nodes that are eventually safe in the graph.
/// \return A vector containing all nodes that are eventually safe.
/// https://leetcode.com/problems/find-eventual-safe-states/
vector<int> Graph::eventualSafeNodes()
{
    unordered_map<int,bool> safe; /// Map to store safe status of nodes.
    vector<int> result; /// Vector to store eventual safe nodes.

    /// Iterate through each node in the graph.
    for (int i = 0; i < this->n; i++)
    {
        /// If the current node is eventually safe, add it to the result vector.
        if (dfsSafeNodes(i, safe))
            result.push_back(i);
    }
    return result; /// Return the vector containing eventual safe nodes.
}

///Public function to find the shortest bridge connecting the two islands in the grid.
/// \param grid Input grid representing the matrix.
/// \return The minimum number of flips needed to connect the two islands.
/// https://leetcode.com/problems/shortest-bridge/
int Graph::shortestBridge(vector<vector<int>>& grid)
{
    int n = grid.size();
    vector<vector<int>> direct = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};
    set<pair<int, int>> visit;

    /// DFS to find the first island and mark its cells as visited
    bool foundFirstIsland = false;
    for (int i = 0; i < n && !foundFirstIsland; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (grid[i][j] == 1)
            {
                dfsBridge(n, i, j, visit, grid, direct);
                foundFirstIsland = true;
                break;
            }
        }
    }

    int res = 0; /// Variable to store the minimum number of flips needed
    bfsBridge(n, visit, grid, direct, res); /// BFS to find the shortest bridge between islands

    return res; /// Return the minimum number of flips required to connect the two islands
}

/// Public function to check if given equations are possible.
/// \param equations Vector of strings representing equations.
/// \return True if the equations are consistent, false otherwise.
/// https://leetcode.com/problems/satisfiability-of-equality-equations/
bool Graph::equationsPossible(vector<string>& equations)
{
    /// Initialize the graph with empty adjacency lists
    this->adjList.resize(128); /// ASCII code of small letters 'a' to 'z'

    /// Build the graph from '==' equations
    for (const string& eq : equations)
    {
        if (eq[1] == '=')
        {
            int x = eq[0];
            int y = eq[3];
            /// Add edges for equality relationships in both directions
            this->adjList[x].push_back(y);
            this->adjList[y].push_back(x);
        }
    }

    /// Check '!=' equations, if they violate the connected components
    for (const string& eq : equations)
    {
        if (eq[1] == '!')
        {
            int x = eq[0];
            int y = eq[3];

            /// Handle special case where the equation is x != x
            if (x == y)
            {
                return false; /// This is always false, as it's a contradiction
            }

            /// Use DFS to check if there's a path between x and y
            vector<bool> visited(128, false);
            if (dfsEquations(x, y, visited))
            {
                return false; /// Conflicting equations found, return false
            }
        }
    }

    return true; /// All equations are consistent, return true
}

/// Public function to find the lowest cost path between the prince and the castle.
/// \return The lowest cost path
/// https://www.infoarena.ro/problema/padure
int Graph::padure(int N, int M, int pl, int pc, int cl, int cc, vector<vector<int>>& forest)
{
    vector<vector<int>> direct = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};
    int result = bfsForest(forest, N, M, {pl, pc}, {cl, cc}, direct);
    return result;
}

/// Public function to determine the orientation of edges in a graph to avoid paths of length two or greater.
/// \param orderList A vector containing the order of edges as they were given in the input.
/// https://codeforces.com/contest/1144/problem/F
void Graph::orientation(vector<int>& orderList)
{
    vector<int> colors(this->n + 1, -1); /// Initialize colors array with -1 (indicating color not assigned).

    /// Check if the graph is bipartite using BFS algorithm.
    if (isBipartite(colors))
    {
        /// If the graph is bipartite, print "YES" to indicate it's possible to orient the edges.
        cout << "YES" << endl;

        /// Iterate through the edges in the order they were given in the input.
        for (int i = 0; i < orderList.size(); i++)
        {
            /// Get the first node of the current edge from the orderList.
            int u = orderList[i];

            /// Determine the orientation of the current edge based on the color of node u.
            /// If u is colored with 0, print "1" (indicating reverse orientation).
            /// If u is colored with 1, print "0" (indicating original orientation).
            if (colors[u] == 0)
            {
                cout << "1";
            }
            else
            {
                cout << "0";
            }
        }
    }
    else
    {
        /// If the graph is not bipartite, print "NO" as it's not possible to orient the edges.
        cout << "NO";
    }
}

/// Public function to find the minimum time required to escape considering all patrol constraints.
/// \param dist: Reference to the 2D vector representing distances between nodes at different times.
/// \return The minimum time required to escape or -1 if not possible.
/// https://www.infoarena.ro/problema/patrol2
int Graph::patrol2(vector<vector<int>>& dist)
{
    /// Perform BFS traversal starting from node 0 and update the distances in 'dist'
    bfsPatrol2(0, dist);

    /// Initialize the minimum distance to a large value
    int minDist = 1e9;

    /// Iterate through all possible times (0 to 419) for the final node
    for (int i = 0; i < 420; i++)
    {
        /// If the final node is reachable at the current time and the time is smaller than the minimum distance
        if (dist[this->n - 1][i] != -1 && dist[this->n - 1][i] < minDist)
        {
            /// Update the minimum distance
            minDist = dist[this->n - 1][i];
        }
    }

    /// If a valid path exists, return the minimum time, otherwise return -1
    if (minDist != 1e9)
        return minDist;
    else
        return -1;
}


/// Public function to calculate the minimum maximum distance among marked nodes in the graph.
/// \param k The number of marked nodes.
/// \param marked A boolean vector indicating whether a node is marked or not.
/// \return The minimum maximum distance among the marked nodes.
/// https://codeforces.com/contest/1881/problem/F
int Graph::minimumMaximum(int k, vector<bool>& marked)
{
    int result; /// Declare a variable to store the minimum maximum distance.

    if (k > 1) /// Check if there are at least two marked nodes for which the minimum maximum distance needs to be calculated.
    {
        vector<int> dist = bfsMinimumMaximum(1); /// Perform BFS from the first node to calculate distances to all nodes.
        int maximum = 0, markedNode;

        for (int i = 1; i <= this->n; i++) /// Find the marked node with the maximum distance from the starting node.
        {
            if (marked[i] && dist[i] > maximum)
            {
                maximum = dist[i];
                markedNode = i;
            }
        }

        maximum = 0; /// Reset the maximum distance.

        dist = bfsMinimumMaximum(markedNode); /// Perform BFS from the marked node with the maximum distance to calculate distances to all nodes.
        for (int i = 1; i <= this->n; i++) /// Find the maximum distance among the marked nodes.
        {
            if (marked[i] && dist[i] > maximum)
                maximum = dist[i];
        }

        result = (maximum + 1) / 2; /// Calculate the minimum maximum distance as half of the maximum distance found.
    }
    else
    {
        result = 0; /// If there is only one marked node, the minimum maximum distance is 0.
    }

    return result; /// Return the calculated minimum maximum distance.
}


int main()
{
    int cnt=0;
    while(true)
    {
        displayMenu();
        int command;
        cin>>command;
        switch(command)
        {
        case 1:
        {
            clearScreen();
            cout<<"n= ";
            int n;
            cin>>n;
            int x, y;
            cout<<endl<<"dislikes (enter -1 to stop)= ";
            vector<vector<int>> dislikes;
            while (cin>>x && x!=-1 && cin>>y)
            {
                dislikes.push_back({x, y});
            }
            cout<<endl;
            Graph G (n,dislikes,false);
            cout<<G.possibleBipartition();
            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;


        }

        case 2:
        {
            clearScreen();
            cout<<"numCourses= ";
            int numCourses;
            cin>>numCourses;
            int x, y;
            cout<<endl<<"prerequisites (enter -1 to stop)= ";
            vector<vector<int>> prerequisites;
            while (cin>>x && x!=-1 && cin>>y)
            {
                prerequisites.push_back({y, x});
            }
            cout<<endl;
            Graph G (numCourses,prerequisites,true);
            vector<int> orderG = G.findOrder();
            for (int i = 0; i < orderG.size(); i++)
            {
                cout << orderG[i] << " ";
            }

            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }

        case 3:
        {
            clearScreen();
            cout<<"n= ";
            int n;
            cin>>n;
            int x, y;
            cout<<endl<<"connections (enter -1 to stop)= ";
            vector<vector<int>> connections;
            while (cin>>x && x!=-1 && cin>>y)
            {
                connections.push_back({x, y});
            }
            cout<<endl;
            Graph G (n,connections,false);
            vector<vector<int>> criticalConnections = G.criticalConnections();
            for (int i = 0; i < criticalConnections.size(); i++)
            {
                for (int j = 0; j < criticalConnections[i].size(); j++)
                    cout<<criticalConnections[i][j]<< " ";
                cout<<endl;
            }
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');

            clearScreen();
            break;
        }
        case 4:
        {
            clearScreen();
            int x, y;
            cout<<"graph (enter -1 to stop)= ";
            vector<vector<int>> graph;
            while (cin>>x && x!=-1 && cin>>y)
            {
                graph.push_back({x, y});
            }
            cout<<endl;
            Graph G (graph.size(),graph,true);
            vector<int> safeNodes = G.eventualSafeNodes();
            for (int i = 0; i < safeNodes.size(); i++)
            {
                cout << safeNodes[i] << " ";
            }

            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }
        case 5:
        {
            clearScreen();
            cout<<"n= ";
            int n;
            cin>>n;
            cout<<endl<<"grid= ";
            vector<vector<int>> grid(n, vector<int>(n, 0));
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    cin >> grid[i][j];
                }
            }
            cout<<endl;
            Graph G;
            cout<<G.shortestBridge(grid);
            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }
        case 6:
        {
            clearScreen();
            cout<<"equations (enter 'stop' to stop)= ";
            vector<string> equations;

            string equation;
            while (true)
            {
                cin >> equation;
                if (equation == "stop")
                {
                    break;
                }
                equations.push_back(equation);
            }
            cout<<endl;
            Graph G;
            cout<<G.equationsPossible(equations);
            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }
        case 7:
        {
            clearScreen();

            int N, M, pl, pc, cl, cc;

            cout << "N= ";
            cin >> N;
            cout << "M= ";
            cin >> M;
            cout << "pl, pc= ";
            cin >> pl >> pc;
            cout << "cl, cc= ";
            cin >> cl >> cc;

            vector<vector<int>> forest(N + 1, vector<int>(M + 1));

            cout << "forest= ";
            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= M; j++)
                {
                    cin >> forest[i][j];
                }
            }
            cout<<endl;
            Graph G;
            cout<<G.padure(N,M,pl,pc,cl,cc,forest);
            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');

            clearScreen();
            break;
        }
        case 8:
        {
            clearScreen();


            int N,M,K;
            cout<<"N,M,K= ";
            cin>>N>>M>>K;


            vector<vector<int>> edges;
            cout<<"edges= ";
            while (M--)
            {
                int x, y;
                cin >> x >> y;
                edges.push_back({x, y});
            }

            for(int i=0; i<N; i++)
                edges.push_back({i,i});

            Graph G(N, edges, false);

            vector<vector<int>> dist(N, vector<int>(420, 1e9));
            for (int i = 0; i < K; i++)
            {
                int L;
                cout<<"L"<<i<<"= ";
                cin >> L;
                for (int j = 0; j < L; j++)
                {
                    int x;
                    cout<<"H"<<i<<j<<"= ";
                    cin >> x;
                    for (int z = j; z < 420; z += L)
                        dist[x][z] = -1;
                }
            }


            cout<<endl;
            cout<<G.patrol2(dist);
            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');

            clearScreen();
            break;
        }
        case 9:
        {
            clearScreen();
            cout<<"n= ";
            int n;
            cin>>n;
            cout<<endl;
            cout<<"m= ";
            int m;
            cin>>m;
            int x, y;
            cout<<endl<<"edges= ";
            vector<vector<int>> edges;
            vector<int> orderList;
            while (m--)
            {
                cin>>x>>y;
                edges.push_back({x, y});
                orderList.push_back(x);
            }
            cout<<endl;
            Graph G (n,edges,false);
            G.orientation(orderList);
            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }
        case 10:
        {
            clearScreen();
            int t;
            cout<<"t= ";
            cin >> t;

            while (t--)
            {
                int n, k;
                cout<<"n, k= ";
                cin >> n >> k;
                int k1 = k;

                vector<bool> marked(n + 1, false);
                cout<<"marked nodes= ";
                while (k1--)
                {
                    int x;
                    cin >> x;
                    marked[x] = true;
                }

                vector<vector<int>> edges;
                cout<<"edges= ";
                int n1 = n - 1;
                while (n1--)
                {
                    int x, y;
                    cin >> x >> y;
                    edges.push_back({x, y});
                }

                Graph G(n, edges, false);
                cout << G.minimumMaximum(k, marked) << endl<<endl;
            }

            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }

        case 11:
        {
            clearScreen();
            cnt=1;
            break;

        }

        }
        if(cnt==1)
        {
            clearScreen();
            break;
        }
    }

    return 0;
}
