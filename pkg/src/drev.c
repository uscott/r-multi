#include "multi.h"

void drev(long len, const double *in, double *out)
{
    out += len - 1;
    while (len-- > 0) *out-- = *in++;
}
