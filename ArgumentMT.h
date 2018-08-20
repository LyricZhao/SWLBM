# define MT_NUMS 64

# define MAXI 3
# define MAXJ 4
# define MAXK 50
# define MAXL 19

# define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
# define BLOCK_HIH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
# define BLOCK_SIZE(id, p, n) (BLOCK_HIH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
# define BLOCK_OWNER(j, p, n) (((p) * ((j) + 1) - 1) / (n))

const Real wM[20] = {
		(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.),0}; // Add an 0

//                     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
const Real e_xM[20] = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0}; // Add an 0
const Real e_yM[20] = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0}; // Add an 0
const Real e_zM[20] = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0}; // Add an 0
