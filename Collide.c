#include "Argument.h"
#include "myArgument.h"
#include <math.h>
/*------------------------------------
 *      Main Computing Part
 *-----------------------------------*/

void collide(Real *****nodes,
		     int ***flags,
			 int Xst,
			 int Xed,
			 int Yst,
			 int Yed,
			 int nz,
			 int current)
{
	int i, j, k, l;
    Real rho, u_x, u_y, u_z;
    Real fi;
	Real feq[19], Qo, omegaNew;
	int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	for(i = Xst; i < Xed; i++) {
		for(j = Yst; j < Yed; j++) {
			for(k = 0; k < nz; k++) {
				if(flags[i - Xst + 1][j - Yst + 1][k] == FLUID || flags[i - Xst + 1][j - Yst + 1][k] == BOUNCE) {
					rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
					for (l = 0; l < 19; l++) {
						fi = nodes[current][i - Xst + 1][j - Yst + 1][k][l];
						rho += fi;
						u_x += e_xM[l] * fi;
						u_y += e_yM[l] * fi;
						u_z += e_zM[l] * fi;
						if(myrank == 0) printf("%.10lf\n", fi);
					}

					u_x /= rho;
					u_y /= rho;
					u_z /= rho;
					if(myrank == 0) printf("%.10f %.10f %.10f %.10f\n", rho, u_x, u_y, u_z);

					// l2
					for (l = 0; l < 19; l++) {
						const Real tmp = (e_xM[l] * u_x + e_yM[l] * u_y + e_zM[l] * u_z);
						feq[l] = wM[l] * rho * (1.0 -
							(3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
							(3.0 *     tmp) +
							(9.0/2.0 * tmp * tmp));
					}
					if(myrank == 0) {
						for(l = 0; l < 19; ++ l) printf("feq: %.20lf\n", feq[l]);
					}

					Qo=0;
					Real Sij[3][3],S;
					Real e[19][3];
					int x1, y1, k1;
					for(k1 = 0; k1 < 19; k1++) {
						e[k1][0] = e_xM[k1];
						e[k1][1] = e_yM[k1];
						e[k1][2] = e_zM[k1];
					}

					for(x1 = 0; x1 < 3; x1++) {
						for(y1 = 0; y1 < 3; y1++) {
							Sij[x1][y1] = 0;
							for(k1 = 0; k1 < 19; k1++) {
								Sij[x1][y1] += e[k1][x1] * e[k1][y1] * (nodes[current][i - Xst + 1][j - Yst + 1][k][k1] - feq[k1]);
							}
							Qo += Sij[x1][y1] * Sij[x1][y1];
						}
					}

					if(myrank == 0) printf("%.10f %.10f %.10f\n", Qo, omega, CSmago);

					nu = (2.0 / omega - 1.0) / 6.0;
					if(myrank == 0) printf("%.10lf\n", nu);
					S = (-nu + sqrt(nu * nu + 18 * CSmago * CSmago * sqrt(Qo))) / 6.0 / CSmago / CSmago;
					if(myrank == 0) printf("%.10lf\n", S);
					omegaNew = 1.0 / (3.0 * (nu + CSmago * CSmago * S) + 0.5);
					if(myrank == 0) printf("%.10f\n", omegaNew);

					for (l = 0; l < 19; l++) {
						nodes[current][i - Xst + 1][j - Yst + 1][k][l] =
							(1.0 - omegaNew) * nodes[current][i - Xst + 1][j - Yst + 1][k][l] +
							omegaNew * feq[l];
					}
					return;
				}
			}
		}
	}
}
