# include <math.h>
# include "athread.h"

typedef float Real;

extern SLAVE_FUN(computeOneStepParallel)();

void computeOneStep(Real *****nodes, char ****walls, int ***flags, int Xst, int Xed, int Yst, int Yed, int nz, int current) {
	long para[9] = {(long) nodes, (long) walls, (long) flags, Xst, Xed, Yst, Yed, nz, current};
	athread_spawn64(SLAVE_FUN(computeOneStepParallel), &para);
	athread_join64();
}
