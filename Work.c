#include "Argument.h"
#include "athread.h"
/*------------------------------------
 *      Main Computing Part
 *-----------------------------------*/
extern SLAVE_FUN(workParallel)();

void work(Real *****nodes,
		    int ****walls,
		    int current)
{
	long t[3]={nodes,walls,current};
	athread_spawn(workParallel,&t);
	athread_join();
}