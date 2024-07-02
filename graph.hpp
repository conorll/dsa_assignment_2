#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <queue>
#include <set>
#include <unordered_map>
#include <limits>
#include <stdexcept>
#include "my_integer.hpp"
#include "index_pq.hpp"

template <typename T>
class Graph {
 private:
  std::vector<std::unordered_map<int, T> > adjList {};
  int numVertices {};

 public:
  // empty graph with N vertices
  explicit Graph(int N);

  // construct graph from edge list in filename
  explicit Graph(const std::string& filename);

  // add an edge directed from vertex i to vertex j with given weight
  void addEdge(int i, int j, T weight);

  // removes edge from vertex i to vertex j
  void removeEdge(int i, int j);

  // is there an edge from vertex i to vertex j?
  bool isEdge(int i, int j) const;

  // return weight of edge from i to j
  // will throw an exception if there is no edge from i to j
  T getEdgeWeight(int i, int j) const;

  // returns number of vertices in the graph
  int size() const;

  // alias a const iterator to our adjacency list type to iterator
  using iterator = 
  typename std::vector<std::unordered_map<int, T> >::const_iterator;

  // cbegin returns const iterator pointing to first element of adjList
  iterator begin() const {
    return adjList.cbegin();
  }

  iterator end() const {
    return adjList.cend();
  }

  // return iterator to a particular vertex
  iterator neighbours(int a) const {
    return adjList.begin() + a;
  }
};

template <typename T>
Graph<T>::Graph(int N) : adjList(N), numVertices {N} {}

template <typename T>
Graph<T>::Graph(const std::string& inputFile) {
  std::ifstream infile {inputFile};
  if (!infile) {
    std::cerr << inputFile << " could not be opened\n";
    return;
  }
  // first line has number of vertices
  infile >> numVertices;
  adjList.resize(numVertices);
  int i {};
  int j {};
  double weight {};
  // assume each remaining line is of form
  // origin dest weight
  while (infile >> i >> j >> weight) {
    addEdge(i, j, static_cast<T>(weight));
  }
}

template <typename T>
int Graph<T>::size() const {
  return numVertices;
}

template <typename T>
void Graph<T>::addEdge(int i, int j, T weight) {
  if (i < 0 or i >= numVertices or j < 0 or j >= numVertices) {
    throw std::out_of_range("invalid vertex number");
  }
  adjList[i].insert({j, weight});
}

template <typename T>
void Graph<T>::removeEdge(int i, int j) {
  // check if i and j are valid
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    adjList[i].erase(j);
  }
}

template <typename T>
bool Graph<T>::isEdge(int i, int j) const {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    return adjList.at(i).contains(j);
  }
  return false;
}

template <typename T>
T Graph<T>::getEdgeWeight(int i, int j) const {
  return adjList.at(i).at(j);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Graph<T>& G) {
  for (int i = 0; i < G.size(); ++i) {
    out << i << ':';
    for (const auto& [neighbour, weight] : *(G.neighbours(i))) {
      out << " (" << i << ", " << neighbour << ")[" << weight << ']';
    }
    out << '\n';
  }
  return out;
}
// return a number of type T to stand in for "infinity"
template <typename T>
T infinity() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  } else if constexpr (std::is_same_v<T, MyInteger>) {
    return MyInteger {std::numeric_limits<int>::max()};
  } else {
    return std::numeric_limits<T>::max();
  }
}

template <typename T>
bool isSubgraph(const Graph<T>& H, const Graph<T>& G) {
  for (int i = 0; i < H.size(); ++i) {

    if (G.neighbours(i) == G.end()) {
      return false;
    }

    auto it = H.neighbours(i);

    for (auto& [neighbor, weight] : *it) {

      if (!G.isEdge(i, neighbor)) {
        return false;
      }

      if (G.getEdgeWeight(i, neighbor) != H.getEdgeWeight(i, neighbor)) {
        return false;
      }
    }
  }

  return true;
}

template <typename T>
bool isTree(const Graph<T>& G, std::vector<bool>& visited, int i) {
  if (visited[i]) {
    return false;
  }

  visited[i] = true;

  auto it = G.neighbours(i);

  for (auto& [neighbor, weight] : *it) {
    if (!isTree(G, visited, neighbor)) {
      return false;
    }
  }

  return true;
}

template <typename T>
bool isTreePlusIsolated(const Graph<T>& G, int root) {
  std::vector<bool> visited(G.size());
  if(!isTree(G, visited, root)) {
    return false;
  }

  for (size_t i = 0; i < visited.size(); ++i) {
    if (visited[i] == false) {
      auto it = G.neighbours(i);

      if ((*it).size() != 0) {
        return false;
      }
    }
  }

  return true;
}

template <typename T>
void dfs_count(const Graph<T>& G, std::vector<T>& bestDistanceTo, int i, T distance) {

  bestDistanceTo[i] = distance;

  auto it = G.neighbours(i);

  for (auto& [neighbor, weight] : *it) {
    dfs_count(G, bestDistanceTo, neighbor, distance + weight);
  }
}

template <typename T>
std::vector<T> pathLengthsFromRoot(const Graph<T>& tree, int root) {
  std::vector<T> bestDistanceTo(tree.size(), infinity<T>());  
  dfs_count(tree, bestDistanceTo, root, T{});
  return bestDistanceTo;
}

template <typename T>
bool allEdgesRelaxed(const std::vector<T>& bestDistanceTo, const Graph<T>& G, int source) {

  if (bestDistanceTo[source] != T{}) {
    return false;
  }

  for (int i = 0; i < G.size(); ++i) {

    auto it = G.neighbours(i);

    for (auto& [neighbor, weight] : *it) {
      if (bestDistanceTo[neighbor] > bestDistanceTo[i] + weight) {
        return false;
      }
    }
  }
  return true;
}

// End of functions from Graph class


// lazy solution as in Tutorial Week 10
template <typename T>
std::vector<T> singleSourceLazyDistance(const Graph<T>& G, int source) {
  // alias the long name for a minimum priority queue holding
  // objects of type DistAndVertex
  using DistAndVertex = std::pair<T, int>;
  using minPQ = std::priority_queue<DistAndVertex,
                                  std::vector<DistAndVertex>,
                                  std::greater<DistAndVertex> >;
  minPQ queue {};
  queue.push({T {}, source});
  // record best distance to vertex found so far
  int N = G.size();
  std::vector<T> bestDistanceTo(N, infinity<T>());
  bestDistanceTo.at(source) = T {};
  // being in visited means we have already explored a vertex's neighbours
  // the bestDistanceTo for a vertex in visited is the true distance.
  std::vector<bool> visited(N);
  while (!queue.empty()) {
    auto [dist, current] = queue.top();
    queue.pop();
    // as we use a lazy version of Dijkstra a vertex can appear multiple
    // times in the queue.  If we have already visited the vertex we
    // take out of the queue we just go on to the next one
    if (visited.at(current)) {
      continue;
    }
    visited.at(current) = true;
    // relax all outgoing edges of current
    for (const auto& [neighbour, weight] : *(G.neighbours(current))) {
      T distanceViaCurrent = bestDistanceTo.at(current) + weight;
      if (bestDistanceTo.at(neighbour) > distanceViaCurrent) {
        bestDistanceTo.at(neighbour) = distanceViaCurrent;
        // lazy dijkstra: nextPoint could already be in the queue
        // we don't update it with better distance just found.
        queue.push({distanceViaCurrent, neighbour});
      }
    }
  }
  return bestDistanceTo;
}

// Implement your lazy solution using std::priority_queue here
template <typename T>
Graph<T> singleSourceLazy(const Graph<T>& G, int source) {
  using DistAndVertex = std::pair<T, int>;
  using minPQ = std::priority_queue<DistAndVertex,
                                  std::vector<DistAndVertex>,
                                  std::greater<DistAndVertex> >;

  int N = G.size();

  std::vector<T> bestDistanceTo(N, infinity<T>());
  std::vector<int> prev(N, -1);
  bestDistanceTo[source] = T{};

  std::vector<bool> visited(N);
  
  minPQ queue {};
  queue.push({T{}, source});

  while (!queue.empty()) {
    const int current = queue.top().second;
    queue.pop();

    if (visited[current]) {
      continue;
    }

    visited[current] = true;

    for (const auto& [neighbour, weight] : *(G.neighbours(current))) {
      T distanceViaCurrent = bestDistanceTo[current] + weight;
      if (distanceViaCurrent < bestDistanceTo[neighbour]) {
        bestDistanceTo[neighbour] = distanceViaCurrent;
        queue.push({distanceViaCurrent, neighbour});
        prev[neighbour] = current;
      }
    }
  }

  Graph<T> shortestPaths(N);

  for (int i = 0; i < N; ++i) {

    int p = prev[i];

    if (p != -1) {
      shortestPaths.addEdge(p, i, G.getEdgeWeight(p, i));
    }
  }


  return shortestPaths;
}


// Implement your solution using an index priority queue here
template <typename T>
Graph<T> singleSourceIndex(const Graph<T>& G, int source) {

  int N = G.size();

  std::vector<bool> visited(N);


  std::vector<T> bestDistanceTo(N, infinity<T>());
  bestDistanceTo[source] = T{};

  std::vector<int> prev(N, -1);

  IndexPriorityQueue<T> queue(N);
  queue.push(T{}, source);

  while (!queue.empty()) {
    const int& current = queue.top().second;

    if (visited[current]) {
      continue;
    }

    visited[current] = true;

    for (const auto& [neighbour, weight] : *(G.neighbours(current))) {
      T distanceViaCurrent = bestDistanceTo[current] + weight;
      if (distanceViaCurrent < bestDistanceTo[neighbour]) {
        bestDistanceTo[neighbour] = distanceViaCurrent;
        queue.changeKey(distanceViaCurrent, neighbour);
        prev[neighbour] = current;
      }
    }
    queue.pop();
  }

  Graph<T> shortestPaths(N);

  for (int i = 0; i < N; ++i) {

    int p = prev[i];

    if (p != -1) {
      shortestPaths.addEdge(p, i, G.getEdgeWeight(p, i));
    }
  }
  return shortestPaths;
}

// Implement your solution using std::set here
template <typename T>
Graph<T> singleSourceSet(const Graph<T>& G, int source) {

  using DistAndVertex = std::pair<T, int>;
  using minPQ = std::set<DistAndVertex, std::greater<DistAndVertex> >;

  int N = G.size();


  std::vector<T> bestDistanceTo(N, infinity<T>());
  bestDistanceTo[source] = T{};

  std::vector<int> prev(N, -1);

  minPQ queue {};
  queue.insert({T{}, source});

  while (!queue.empty()) {

    const int current = queue.begin()->second;
    queue.erase(queue.begin());

    for (const auto& [neighbour, weight] : *(G.neighbours(current))) {
      T distanceViaCurrent = bestDistanceTo[current] + weight;
      if (distanceViaCurrent < bestDistanceTo[neighbour]) {
        bestDistanceTo[neighbour] = distanceViaCurrent;
        queue.insert({distanceViaCurrent, neighbour});
        prev[neighbour] = current;
      }
    }
  }
  Graph<T> shortestPaths(N);

  for (int i = 0; i < N; ++i) {

    int p = prev[i];

    if (p != -1) {
      shortestPaths.addEdge(p, i, G.getEdgeWeight(p, i));
    }
  }
  return shortestPaths;

}

// put your "best" solution here
// this is the one we will use for performance testing
template <typename T>
Graph<T> singleSourceShortestPaths(const Graph<T>& G, int source) {
  return singleSourceIndex(G, source);
}

#endif      // GRAPH_HPP_




