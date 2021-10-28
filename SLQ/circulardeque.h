#ifndef CIRCULARDEQUE_H
#define CIRCULARDEQUE_H
#include <cassert>
#include <cstddef>
#include <utility>

template <typename T>
class CircularDeque {
  size_t n, l, r;
  bool isFull;
  T *arr;

public:
  CircularDeque(size_t n) : n{n}, l{0}, r{0}, isFull{0} {
    arr = new T[n];
  }

  CircularDeque &operator=(const CircularDeque<T> &Q) {
    n = Q.n; l = Q.l; r = Q.r;
    isFull = Q.isFull;
    delete[] arr;
    arr = Q.arr;
    return *this;
  }

  CircularDeque &operator=(CircularDeque<T> &&Q) {
    n = Q.n; l = Q.l; r = Q.r;
    isFull = Q.isFull;
    delete[] arr;
    arr = Q.arr;
    Q.arr = nullptr;
    return *this;
  }

  ~CircularDeque() { delete[] arr; }

  bool isempty() const { return !isFull && l == r; }

  void clear() { l = r = isFull = 0; }

  size_t capacity() const { return n; }

  size_t length() const { return isFull ? n : (r >= l ? r - l : r + n - l); }

  void push(const T &t) {
    assert(!isFull);
    arr[r] = t;
    r = (r == n - 1 ? 0 : r + 1);
    if (l == r) {
      isFull = 1;
    }
  }

  void push(T &&t) {
    assert(!isFull);
    arr[r] = t;
    r = (r == n - 1 ? 0 : r + 1);
    if (l == r) {
      isFull = 1;
    }
  }

  T pop() {
    if (l == r) {
      isFull = 0;
    }
    r = (r == 0 ? n - 1 : r - 1);
    return arr[r];
  }

  void pushfirst(const T &t) {
    assert(!isFull);
    l = (l == 0 ? n - 1 : l - 1);
    if (l == r) {
      isFull = 1;
    }
    arr[l] = t;
  }

  void pushfirst(T &&t) {
    assert(!isFull);
    l = (l == 0 ? n - 1 : l - 1);
    if (l == r) {
      isFull = 1;
    }
    arr[l] = t;
  }

  T popfirst() {
    if (l == r) {
      isFull = 0;
    }
    size_t prev = l;
    l = (l == n - 1 ? 0 : l + 1);
    return arr[prev];
  }

  T &first() { return arr[l]; }

  const T &first() const { return arr[l]; }

  T &last() { return arr[r == 0 ? n - 1 : r - 1]; }

  const T &last() const { return arr[r == 0 ? n - 1 : r - 1]; }
};
#endif // CIRCULARDEQUE_H
