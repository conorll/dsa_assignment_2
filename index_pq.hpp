#ifndef INDEX_PRIORITY_QUEUE_HPP_
#define INDEX_PRIORITY_QUEUE_HPP_

#include <vector>
#include <algorithm>
#include <iostream>


template <typename T>
class IndexPriorityQueue {
 private:
  // vector to hold priorities.  
  // priorities.at(i) is the priority associated to index i
  // the heap will only contain indices, we look up priorities
  // in this vector as needed
  std::vector<T> priorities {};
  // priorityQueue stores indices: priorityQueue.at(i) is an index
  // priorityQueue functions as the heap and is heap ordered: 
  // priorities.at(priorityQueue.at(i)) <= priorities.at(priorityQueue.at(2 * i)) 
  // priorities.at(priorityQueue.at(i)) <= priorities.at(priorityQueue.at(2 * i) + 1) 
  std::vector<int> priorityQueue {};
  // indexToPosition.at(i) is the position in priorityQueue of index i
  // priorityQueue.at(indexToPosition.at(i)) = i
  // indexToPosition.at(priorityQueue.at(j)) = j
  std::vector<int> indexToPosition {};
  int size_ = 0;

 public:
  explicit IndexPriorityQueue(int);
  void push(const T&, int);
  void pop();
  void erase(int);
  bool contains(int) const;
  void changeKey(const T&, int);
  std::pair<T, int> top() const;
  bool empty() const;
  int size() const;
  void output();

 private:
  // TODO: you may want to add your own member functions. swim? sink?
  void swim(int i);
  void sink(int i);
  void swap(int i, int j);
};

// Useful helper functions
int leftChild(int i) {
  return 2*i + 1;
}

int rightChild(int i) {
  return 2*i + 2;
}

int parent(int i) {

  if (i == 0) {
    return 0;
  }

  return (i - 1)/2;
}

template <typename T>
void IndexPriorityQueue<T>::output() {
  for (int i = 0; i < priorities.size(); ++i) {
    std::cout << priorities[i];
  }
  std::cout << "\n";
  for (int i = 0; i < priorities.size(); ++i) {
    std::cout << indexToPosition[i];
  }
  std::cout << "\n";
  for (int i = 0; i < priorities.size(); ++i) {
    std::cout << priorityQueue[i];
  }
  std::cout << "\n";
  std::cout << "\n";
}

// IndexPriorityQueue member functions
template <typename T>
IndexPriorityQueue<T>::IndexPriorityQueue(int N) {
  priorities.resize(N, T{});
  indexToPosition.resize(N, -1);
  priorityQueue.resize(N, -1);
}

template <typename T>
bool IndexPriorityQueue<T>::empty() const {
  return size_ == 0;
}

template <typename T>
int IndexPriorityQueue<T>::size() const {
  return size_;
}

template <typename T>
void IndexPriorityQueue<T>::push(const T& priority, int i) {
  if (i < 0 || i >= static_cast<int>(indexToPosition.size()) || indexToPosition[i] != -1) {
    return;
  }

  priorities[i] = priority;
  indexToPosition[i] = size_;
  priorityQueue[size_] = i;

  ++size_;
  swim(size_ - 1);
}

template <typename T>
void IndexPriorityQueue<T>::swim(int pos) {
  while(priorities[priorityQueue[pos]] < priorities[priorityQueue[parent(pos)]]) {
    swap(pos, parent(pos));
    pos = parent(pos);
  }
}

template <typename T>
void IndexPriorityQueue<T>::pop() {
  
  if (size_ == 0) {
    return;
  }

  erase(priorityQueue[0]);
}

template <typename T>
void IndexPriorityQueue<T>::sink(int pos) {
  
  while (true) {

    int minP = pos;
    if (leftChild(pos) < size_ && priorities[priorityQueue[leftChild(pos)]] < priorities[priorityQueue[minP]]) {
      minP = leftChild(pos);
    }
    if (rightChild(pos) < size_ && priorities[priorityQueue[rightChild(pos)]] < priorities[priorityQueue[minP]]) {
      minP = rightChild(pos);
    }
    if (pos == minP) {
      return;
    }
    
    swap(pos, minP);

    pos = minP;
  }
}

template <typename T>
void IndexPriorityQueue<T>::swap(int pos1, int pos2) {
  std::swap(priorityQueue[pos1], priorityQueue[pos2]);
  std::swap(indexToPosition[priorityQueue[pos1]], indexToPosition[priorityQueue[pos2]]);
}

template <typename T>
void IndexPriorityQueue<T>::erase(int i) {

  if (i >= static_cast<int>(indexToPosition.size()) || i < 0) {
    return;
  }

  int pos = indexToPosition[i];

  if (pos == -1) {
    return;
  }

  swap(pos, size_ - 1);
  indexToPosition[i] = -1;
  priorityQueue[size_ - 1] = -1;
  priorities[i] = T{};
  --size_;

  if (pos < size_) {
    swim(pos);
    sink(pos);
  }
}

template <typename T>
std::pair<T, int> IndexPriorityQueue<T>::top() const {
  return {priorities[priorityQueue[0]], priorityQueue[0]};
}

// if vertex i is not present, insert it with key
// otherwise change the associated key value of i to key
template <typename T>
void IndexPriorityQueue<T>::changeKey(const T& key, int i) {

  if (indexToPosition[i] != -1) {
    priorities[i] = key;

    sink(indexToPosition[i]);
    swim(indexToPosition[i]);
  }
  else {
    push(key, i);
  }
}

template <typename T>
bool IndexPriorityQueue<T>::contains(int i) const {
  
  if (i < 0 || i >= size_) {
    return false;
  }

  return indexToPosition[i] != -1;
}


#endif      // INDEX_PRIORITY_QUEUE_HPP_

