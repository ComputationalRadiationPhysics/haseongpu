#include<assert.h>
#include "dynintset.h"

__device__ void DynIntSet::resize(const size_t newSize) {
  int* newBuffer = new int[newSize];
  bufferSize = newSize;
  fillSize = min((unsigned) fillSize, (unsigned) newSize);

  for(int i=0; i<fillSize; ++i) {
    newBuffer[i] = buffer[i];
  }
  delete[] buffer;
  buffer = newBuffer;
}

__device__ DynIntSet::DynIntSet(int initialSize) {
  buffer = new int[initialSize];
  bufferSize = initialSize;
  fillSize = 0;
}

__device__ DynIntSet::~DynIntSet() {
  delete[] buffer;
}

//int& DynIntSet::operator[](const unsigned i) {
//  assert(i<fillSize);
//  return buffer[i];
//}
__device__ int DynIntSet::operator[](const unsigned i) const {
  assert(i<fillSize);
  return buffer[i];
}

__device__ int DynIntSet::getFillSize() const {
  return fillSize;
}
__device__ int DynIntSet::getBufferSize() const {
  return bufferSize;
}

__device__ void DynIntSet::push(int a) {
  buffer[fillSize++] = a;
}
__device__ void DynIntSet::push(const int *a, const size_t size) {
  if (fillSize + size > bufferSize) resize(bufferSize*2);

  for(int i=0; i<size; ++i) {
    bool present = false;
    for(int j=0; j<fillSize; ++j) {
      if(a[i] == buffer[j]){
        present = true;
        break;
      }
    }
    if(!present) push(a[i]);
  }
}
