#include <cstring>

#include "testing_tool_rq.h"

using namespace std;
using namespace uni10;

int main(int argc, char **argv)
{
  if (argc == 2 && !strcmp(argv[1], "--inplace"))
    testing_rq<uni10_complex128>(true);
  else
    testing_rq<uni10_complex128>(false);
  return 0;
}
