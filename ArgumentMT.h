# include "slave.h"

# define MT_NUMS 64

# define MAXI 3
# define MAXJ 4
# define MAXK 50
# define MAXL 19

# define FLUID 0
# define NOSLIP 1
# define VELOCITY 2
# define BOUNCE 3
# define PRESSURE 4

# define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
# define BLOCK_HIH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
# define BLOCK_SIZE(id, p, n) (BLOCK_HIH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
# define BLOCK_OWNER(j, p, n) (((p) * ((j) + 1) - 1) / (n))

typedef float Real;

__thread_local_fix const Real wM[20] __attribute__((aligned(32))) = {
		(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.), (0.0)}; // Add an 0

__thread_local_fix const Real e_xM[20] __attribute__((aligned(32))) = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0}; // Add an 0
__thread_local_fix const Real e_yM[20] __attribute__((aligned(32))) = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0}; // Add an 0
__thread_local_fix const Real e_zM[20] __attribute__((aligned(32))) = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0}; // Add an 0

__thread_local_fix const int dfInvM[20] __attribute__((aligned(32))) = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18, 0};

__thread_local_fix const int e_xI[20] __attribute__((aligned(32))) = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0}; // Add an 0
__thread_local_fix const int e_yI[20] __attribute__((aligned(32))) = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0}; // Add an 0
__thread_local_fix const int e_zI[20] __attribute__((aligned(32))) = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0}; // Add an 0
