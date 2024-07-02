#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include "graph.hpp"
#include "my_integer.hpp"

// You can generate some random graphs to help in your testing
// The graph has N vertices and p is the probability there is an
// edge between any two vertices. 
// You can vary seed to get different graphs
Graph<int> randomGraph(int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  // set up random number generator that is 1 with probability p and
  // 0 with probability 1-p
  std::binomial_distribution<int> heads {1, p};
  // Set the minimum and maximum edge weight here
  const int minWeight {1};
  const int maxWeight {10};
  std::uniform_int_distribution<int> weight {minWeight, maxWeight};
  Graph<int> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        // add edge between i and j with random weight 
        // between minWeight and maxWeight
        G.addEdge(i, j, weight(mt));
      }
    }
  }
  return G;
}

// Example Test of shortest path algorithm from tutorial
TEST(DijkstraTest, distanceFrom0LazyInt) {
  Graph<int> G {"tinyEWD.txt"};
  std::vector<int> result = singleSourceLazyDistance<int>(G, 0);
  EXPECT_EQ(result[0], 0);
  EXPECT_EQ(result[1], 105);
  EXPECT_EQ(result[2], 26);
  EXPECT_EQ(result[3], 99);
  EXPECT_EQ(result[4], 38);
  EXPECT_EQ(result[5], 73);
  EXPECT_EQ(result[6], 151);
  EXPECT_EQ(result[7], 60);
}

// Example showing how to use MyInteger class to get performance data
TEST(DijkstraTest, distanceFrom0LazyMyInteger) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  MyInteger::clearCounts();
  std::vector<MyInteger> result = singleSourceLazyDistance<MyInteger>(G, 0);
  // see how many operations on MyInteger the algorithm made
  MyInteger::printCounts();
  EXPECT_EQ(result[0], MyInteger {0});
  EXPECT_EQ(result[1], MyInteger {105});
  EXPECT_EQ(result[2], MyInteger {26});
  EXPECT_EQ(result[3], MyInteger {99});
  EXPECT_EQ(result[4], MyInteger {38});
  EXPECT_EQ(result[5], MyInteger {73});
  EXPECT_EQ(result[6], MyInteger {151});
  EXPECT_EQ(result[7], MyInteger {60});
}

// Lazy tiny

TEST(DijkstraTest, singleSourceLazyTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceLazy<MyInteger>(G, 0);
  // see how many operations on MyInteger the algorithm made
  MyInteger::printCounts();
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<MyInteger> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceLazyTinyInt) {
  Graph<int> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<int> result = singleSourceLazy<int>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<int> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceLazyTinyDouble) {
  Graph<double> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<double> result = singleSourceLazy<double>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<double> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

// Lazy medium

TEST(DijkstraTest, singleSourceLazyMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceLazy<MyInteger>(G, 0);
  // see how many operations on MyInteger the algorithm made
  MyInteger::printCounts();
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<MyInteger> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceLazyMediumInt) {
  Graph<int> G {"mediumEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<int> result = singleSourceLazy<int>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<int> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceLazyMediumDouble) {
  Graph<double> G {"mediumEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<double> result = singleSourceLazy<double>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<double> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

// Lazy random

TEST(DijkstraTest, singleSourceLazyRandom) {
  Graph<int> G = randomGraph(400, 43);
  // set all counts to 0 before running shortest path algorithm
  Graph<int> result = singleSourceLazy<int>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<int> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

// Index tiny

TEST(DijkstraTest, singleSourceIndexTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 0);
  // see how many operations on MyInteger the algorithm made
  MyInteger::printCounts();
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<MyInteger> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceIndexTinyInt) {
  Graph<int> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<int> result = singleSourceIndex<int>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<int> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceIndexTinyDouble) {
  Graph<double> G {"tinyEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<double> result = singleSourceIndex<double>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<double> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

// Index medium

TEST(DijkstraTest, singleSourceIndexMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 0);
  // see how many operations on MyInteger the algorithm made
  MyInteger::printCounts();
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<MyInteger> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceIndexMediumInt) {
  Graph<int> G {"mediumEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<int> result = singleSourceIndex<int>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<int> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

TEST(DijkstraTest, singleSourceIndexMediumDouble) {
  Graph<double> G {"mediumEWD.txt"};
  // set all counts to 0 before running shortest path algorithm
  Graph<double> result = singleSourceIndex<double>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<double> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

// Index random

TEST(DijkstraTest, singleSourceIndexRandom) {
  Graph<int> G = randomGraph(400, 43);
  // set all counts to 0 before running shortest path algorithm
  Graph<int> result = singleSourceIndex<int>(G, 0);
  // see how many operations on MyInteger the algorithm made
  EXPECT_EQ(isTreePlusIsolated(result, 0), true);
  EXPECT_EQ(isSubgraph(result, G), true);

  std::vector<int> bestDistanceTo = pathLengthsFromRoot(result, 0);

  EXPECT_EQ(allEdgesRelaxed(bestDistanceTo, G, 0), true);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

