# include "simd.h"

# include <stdio.h>

int main() {

  printf("%u %u %u\n", sizeof(floatv4), sizeof(float), sizeof(int));

  floatv4 sV = 1.0;
  float *f = (float*) &sV;

  int i;
  for(i = 0; i < 4; ++ i) printf("%.2lf ", f[i]); printf("\n");

  return 0;
}
