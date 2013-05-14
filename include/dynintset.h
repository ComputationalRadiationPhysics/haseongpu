#ifndef DYNINTSET_H
#define DYNINTSET_H 

#include <cstddef>

class DynIntSet {
  int *buffer;
  size_t fillSize;
  size_t bufferSize;

  __device__ void resize(const size_t newSize);

  public:
  __device__ DynIntSet(int initialSize);
  __device__ ~DynIntSet();

  //int& operator[](const size_t i);
  __device__ int operator[](const unsigned i) const;

  __device__ int getFillSize() const;
  __device__ int getBufferSize() const;

  __device__ void push(int a);
  __device__ void push(const int *a, const size_t size);
};

#endif /* DYNINTSET_H */
