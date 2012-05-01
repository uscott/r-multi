#include "uri.h"
void dmult(long len, const double *in1, const double *in2, double *out)
{
  while (len-- > 0)
    *out++ = (*in1++) * (*in2++);
}
