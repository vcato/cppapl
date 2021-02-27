#include "optional.hpp"


int main()
{
  assert(Optional<int>(4) != Optional<int>(5));
  assert(Optional<int>() != Optional<int>(4));
  assert(!(Optional<int>() != Optional<int>()));
  assert(Optional<int>(4) == Optional<int>(4));
}
