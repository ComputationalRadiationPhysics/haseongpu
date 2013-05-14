#include "dynintset.h"
#include<assert.h>

int main() {
  DynIntSet set(16);

  for(int i=0; i<15; ++i) set.push(i);
  assert(set.getFillSize() == 15);
  assert(set.getBufferSize() == 16);
  assert(set[0] == 0);
  assert(set[6] == 6);
  assert(set[14] == 14);

  int test2[] = {15,16,17};
  set.push(test2, 3);
  assert(set.getFillSize() == 18);
  assert(set.getBufferSize() == 32);
  assert(set[15] == 15);
  assert(set[17] == 17);

  int test1[] = {1,1,1};
  set.push(test1, 3);
  assert(set.getFillSize() == 18);
}
