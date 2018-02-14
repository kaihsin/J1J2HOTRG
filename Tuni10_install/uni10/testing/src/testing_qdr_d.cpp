#include <cstring>

#include "testing_tool_qdr.h"

using namespace std;
using namespace uni10;

int main(int argc, char **argv)
{
  if (argc == 2 && !strcmp(argv[1], "--inplace"))
    testing_qdr<uni10_double64>(true);
  else
    testing_qdr<uni10_double64>(false);
  return 0;
}
