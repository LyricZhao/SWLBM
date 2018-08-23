#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Argument.h"
#include "athread.h"

extern SLAVE_FUN(computeOneStepParallel)();
extern SLAVE_FUN(preworkSlave)();

inline Real sqr(Real s) {
  return s * s;
}

int main(int argc, char *argv[])
{
    int ***flags;
    int i, j, k, s, l;
    int myrank, my2drank, size;
    int current = 0, other = 1;
    int periods[NUM_DIMS] = {0, 0};
    int dims[NUM_DIMS] = {0, 0};
    int coords[2];
    Real df;
    unsigned int X_section ;
    unsigned int X_res     ;
    unsigned int Y_section ;
    unsigned int Y_res     ;
    int	Xst ;
    int	Xed ;
    int	Yst ;
    int	Yed ;
    int x_sec ;
    int y_sec ;
    Real *****nodes;
    int ****walls;
    Real ***temp_right;
    Real ***temp_left;
    Real ***temp_down;
    Real ***temp_up;
    Real ***temp_right_send;
    Real ***temp_left_send;
    Real ***temp_down_send;
    Real ***temp_up_send;
    Real **temp_lu;
    Real **temp_ld;
    Real **temp_ru;
    Real **temp_rd;
    Real **temp_lu_send;
    Real **temp_ld_send;
    Real **temp_ru_send;
    Real **temp_rd_send;
    MPI_Comm mycomm;
    MPI_Status sta[16];
    MPI_Request req[16];
    int n = 0;
    int count;
    int local_rankinfo[7];
    int **rankinfo;
    Real **image;
    Real **local_image;

	/*------------------------
         * Parameter Set
         * ----------------------*/

    	setParameter();

	/*---------------------------*
	 * MPI Init
	 * --------------------------*/


    	MPI_Init(&argc, &argv);

    	MPI_Comm_size(MPI_COMM_WORLD, &size);

    	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	MLOG("CPC -- LBM-simulation Benchmark!\n\n");

	/*----------------------------*
         * MPI Mesh Section
         * ---------------------------*/

	MPI_Dims_create(size, NUM_DIMS, dims);

    	MPI_Cart_create(MPI_COMM_WORLD, NUM_DIMS, dims, periods, 0, &mycomm);

    	MPI_Comm_rank(mycomm, &my2drank);

    	MPI_Cart_coords(mycomm, my2drank, 2, coords);

	SetMPI(mycomm, dims, coords);



    	X_section = X / dims[0];
    	X_res     = X % dims[0];
    	Y_section = Y / dims[1];
    	Y_res     = Y % dims[1];

	Xst = (coords[0] < X_res) ? (coords[0] * (X_section + 1))
								  : (coords[0] *  X_section + X_res);
	Xed = (coords[0] < X_res) ? (Xst       + (X_section + 1))
								  : (Xst       +  X_section);
	Yst = (coords[1] < Y_res) ? (coords[1] * (Y_section + 1))
								  : (coords[1] *  Y_section + Y_res);
	Yed = (coords[1] < Y_res) ? (Yst       + (Y_section + 1))
								  : (Yst       +  Y_section);

    	x_sec = Xed - Xst;
    	y_sec = Yed - Yst;


	/*-----------------------------*
         *-----------------------------*
         * ----------------------------*/

	MLOG("Size                 :  %d x %d x %d\n", X, Y, Z);
        MLOG("Steps                :  %d\n", STEPS);
        MLOG("Number of Process    :  %d\n\n", size);

        /*-----------------------------*
         * Space Allocate
         * ----------------------------*/



	flags = array3DI(x_sec + 2, y_sec + 2, Z);
	nodes = array5DF(2, x_sec + 2, y_sec + 2, Z, 19);
	walls = array4DI(x_sec, y_sec, Z, 19);
    	temp_lu_send    = array2DF(Z, 19);
    	temp_lu         = array2DF(Z, 19);
    	temp_ld         = array2DF(Z, 19);
    	temp_ru         = array2DF(Z, 19);
    	temp_rd         = array2DF(Z, 19);
    	temp_ld_send    = array2DF(Z, 19);
    	temp_ru_send    = array2DF(Z, 19);
    	temp_rd_send    = array2DF(Z, 19);
    	temp_right      = array3DF(y_sec, Z, 19);
    	temp_left       = array3DF(y_sec, Z, 19);
    	temp_down       = array3DF(x_sec, Z, 19);
    	temp_up         = array3DF(x_sec, Z, 19);
    	temp_right_send = array3DF(y_sec, Z, 19);
    	temp_left_send  = array3DF(y_sec, Z, 19);
    	temp_down_send  = array3DF(x_sec, Z, 19);
    	temp_up_send    = array3DF(x_sec, Z, 19);
    	memset(&walls[0][0][0][0], 0, x_sec * y_sec * Z * 19 * sizeof(int));

	/*----------------------------------------*
	 * local rank info
	 * ---------------------------------------*/

	local_image = array2DF(y_sec, x_sec * 5);

	if(my2drank == 0) {
		image = array2DF(Y, X * 5);
		rankinfo = array2DI(size, 7);
	} else {
		image = array2DF(1, 1);
		rankinfo = array2DI(1, 1);
	}

	INITINPUT(X, Y, Z, Xst, Xed, Yst, Yed, x_sec, y_sec, myrank, size, argv[1], local_rankinfo, rankinfo, flags);

	MLOG("Init >> init Pointer!\n");

    	init_Pointer(flags, nodes, walls, Xst, Xed, Yst, Yed, Z, x_sec, y_sec, 1, LDC_VELOCITY);

	MLOG("Step >> Main Steps start!\n\n");
	/*----------------------------------------------------*
	 * Main Calculation section
	 * ---------------------------------------------------*/
	TIME_ST();
  	athread_init();
  	athread_enter64();
  	char ****wallsChar = (char ****) array4DI(x_sec + 2, y_sec + 2, Z, 5);
	long para[3]={(long)walls,(long)flags,(long)wallsChar};
	athread_spawn64(SLAVE_FUN(preworkSlave), &para);
	athread_join64();

	volatile int updateflag;
	Real rho, u_x, u_y, u_z,fi,feq[19],nfSub[19],Qo, omegaNew,S;
	Real _nu= (2.0 / omega - 1.0) / 6.0;
	for (s = 0; s < STEPS; s++) {
		updateflag=0;
		long para[10] = {(long) nodes, (long) wallsChar, (long) flags, Xst, Xed, Yst, Yed, Z, current, (long) &updateflag};
		athread_spawn64(SLAVE_FUN(computeOneStepParallel), &para);
        bounce_send_init(X,
			 Y,
			 Z,
			 Xst,
		         Xed,
			 Yst,
			 Yed,
			 x_sec,
			 y_sec,
			 current,
			 other,
			 nodes,
                         temp_left_send,
			 temp_right_send,
			 temp_up_send,
			 temp_down_send,
                         temp_ld_send,
			 temp_lu_send,
			 temp_rd_send,
			 temp_ru_send);

        bounce_communicate(mycomm,
		           dims,
			   coords,
			   x_sec,
			   y_sec,
			   Z,
			   &count,
			   sta,
			   req,
			   temp_left_send,
			   temp_right_send,
			   temp_up_send,
			   temp_down_send,
			   temp_left,
			   temp_right,
			   temp_up,
			   temp_down,
			   temp_lu_send,
			   temp_ld_send,
			   temp_ru_send,
			   temp_rd_send,
			   temp_lu,
			   temp_ld,
			   temp_ru,
			   temp_rd);

	for(i = 0; i < count; i++) {
		MPI_Wait(&req[i], &sta[i]);
	}

        bounce_update(X,
		      Y,
		      Z,
		      Xst,
		      Xed,
		      Yst,
		      Yed,
		      myrank,
		      x_sec,
		      y_sec,
		      other,
		      nodes,
                      temp_left,
		      temp_right,
		      temp_up,
		      temp_down,
                      temp_ld,
		      temp_lu,
		      temp_rd,
		      temp_ru);
			updateflag=1;
		{
			for(i = Xst; i < Xed; i++) {
				for(j = Yst; j < Yed; j++) {
					if(i!=Xst&&i!=Xed-1) continue;
					for(k = 0; k < Z; k++) {
						Real *npc0=nodes[current][i - Xst + 1][j - Yst + 1][k];
						Real ****nodesOtherLocal=nodes[other];
						int is1 = i-Xst, in0 = i-Xst+1, ia1 = i-Xst+2;
						int js1 = j-Yst, jn0 = j-Yst+1, ja1 = j-Yst+2;
						if(flags[in0][jn0][k] == 0) {
							npc0[ 0] = nodesOtherLocal[in0][js1][k   ][ 0];
							npc0[ 1] = nodesOtherLocal[in0][ja1][k   ][ 1];
							npc0[ 2] = nodesOtherLocal[is1][jn0][k   ][ 2];
							npc0[ 3] = nodesOtherLocal[ia1][jn0][k   ][ 3];
							npc0[ 4] = nodesOtherLocal[in0][jn0][k - 1][ 4];
							npc0[ 5] = nodesOtherLocal[in0][jn0][k + 1][ 5];
							npc0[ 6] = nodesOtherLocal[is1][js1][k    ][ 6];
							npc0[ 7] = nodesOtherLocal[ia1][js1][k    ][ 7];
							npc0[ 8] = nodesOtherLocal[is1][ja1][k    ][ 8];
							npc0[ 9] = nodesOtherLocal[ia1][ja1][k    ][ 9];
							npc0[10] = nodesOtherLocal[in0][js1][k - 1][10];
							npc0[11] = nodesOtherLocal[in0][js1][k + 1][11];
							npc0[12] = nodesOtherLocal[in0][ja1][k - 1][12];
							npc0[13] = nodesOtherLocal[in0][ja1][k + 1][13];
							npc0[14] = nodesOtherLocal[is1][jn0][k - 1][14];
							npc0[15] = nodesOtherLocal[is1][jn0][k + 1][15];
							npc0[16] = nodesOtherLocal[ia1][jn0][k - 1][16];
							npc0[17] = nodesOtherLocal[ia1][jn0][k + 1][17];
							npc0[18] = nodesOtherLocal[in0][jn0][k    ][18];
						}
						else if(flags[in0][jn0][k] == 3) {
							for(l = 0; l < 19; l++) {
								int inv = dfInv[l];
								if(walls[is1][js1][k][l]) {
									npc0[l] = nodes[other][in0][jn0][k][inv];
								} else {
									npc0[l] = nodes[other][in0 + e_x[inv]][jn0 + e_y[inv]][k + e_z[inv]][l];
								}
							}
						}
						if(flags[in0][jn0][k] == 0 || flags[in0][jn0][k] == 3) {
							rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
							for (l = 0; l < 19; l++) {
								fi = npc0[l];
								rho += fi;
								u_x += e_x[l] * fi;
								u_y += e_y[l] * fi;
								u_z += e_z[l] * fi;
							}

							u_x /= rho;
							u_y /= rho;
							u_z /= rho;

							for (l = 0; l < 19; l++) {
								const Real tmp = (e_x[l] * u_x + e_y[l] * u_y + e_z[l] * u_z);
								feq[l] = w[l] * rho * (1.0 -
									(1.5 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
									(3.0 *     tmp) +
									(4.5 * tmp * tmp));
								nfSub[l]=npc0[l]-feq[l];
							}
							Qo=0;
							Qo += sqr(nfSub[2] + nfSub[3] + nfSub[6] + nfSub[7] + nfSub[8] + nfSub[9] + nfSub[14] + nfSub[15] + nfSub[16] + nfSub[17]);
							Qo += sqr(nfSub[6] - nfSub[7] - nfSub[8] + nfSub[9]);
							Qo += sqr(nfSub[14] - nfSub[15] - nfSub[16] + nfSub[17]);
							double sij = nfSub[6] - nfSub[7] - nfSub[8]+nfSub[9];
							Qo += sij * sij;
							Qo += sqr(nfSub[0] + nfSub[1] + nfSub[6] + nfSub[7] + nfSub[8] + nfSub[9] + nfSub[10] + nfSub[11] + nfSub[12] + nfSub[13]);
							Qo += sqr(nfSub[10] - nfSub[11] - nfSub[12] + nfSub[13]);
							Qo += sqr(nfSub[14] - nfSub[15] - nfSub[16] + nfSub[17]);
							Qo += sqr(nfSub[10] - nfSub[11] - nfSub[12] + nfSub[13]);
							Qo += sqr(nfSub[4] + nfSub[5] + nfSub[10] + nfSub[11] + nfSub[12] + nfSub[13] + nfSub[14] + nfSub[15] + nfSub[16] + nfSub[17]);
							S = (-_nu + sqrt(_nu * _nu + 18 * CSmago * CSmago * sqrt(Qo))) / 6.0 / CSmago / CSmago;
							omegaNew = 1.0 / (3.0 * (_nu + CSmago * CSmago * S) + 0.5);

							rho = 1.0 - omegaNew;
							npc0[ 0] = rho * npc0[ 0] + omegaNew * feq[ 0];
							npc0[ 1] = rho * npc0[ 1] + omegaNew * feq[ 1];
							npc0[ 2] = rho * npc0[ 2] + omegaNew * feq[ 2];
							npc0[ 3] = rho * npc0[ 3] + omegaNew * feq[ 3];
							npc0[ 4] = rho * npc0[ 4] + omegaNew * feq[ 4];
							npc0[ 5] = rho * npc0[ 5] + omegaNew * feq[ 5];
							npc0[ 6] = rho * npc0[ 6] + omegaNew * feq[ 6];
							npc0[ 7] = rho * npc0[ 7] + omegaNew * feq[ 7];
							npc0[ 8] = rho * npc0[ 8] + omegaNew * feq[ 8];
							npc0[ 9] = rho * npc0[ 9] + omegaNew * feq[ 9];
							npc0[10] = rho * npc0[10] + omegaNew * feq[10];
							npc0[11] = rho * npc0[11] + omegaNew * feq[11];
							npc0[12] = rho * npc0[12] + omegaNew * feq[12];
							npc0[13] = rho * npc0[13] + omegaNew * feq[13];
							npc0[14] = rho * npc0[14] + omegaNew * feq[14];
							npc0[15] = rho * npc0[15] + omegaNew * feq[15];
							npc0[16] = rho * npc0[16] + omegaNew * feq[16];
							npc0[17] = rho * npc0[17] + omegaNew * feq[17];
							npc0[18] = rho * npc0[18] + omegaNew * feq[18];
						}
					}
				}
			}
		}
  		//computeOneStep(nodes, wallsChar, flags, Xst, Xed, Yst, Yed, Z, current);
		athread_join64();
		other = current;
		current = current^1;

		if(myrank == 0 && STEPS >= 10 && (s + 1)%(STEPS/10) == 0.0) {
				n += 1;
				MLOG("Step >> [%d/%d] Calculation Completed %d%% \n", s + 1, STEPS, n * 10);
		}

	}

  athread_leave64();
  athread_halt();
  arrayFree4DI((int ****)wallsChar);
	TIME_ED();
	/*-----------------------------*
 	 * OUTPUT
 	 *-----------------------------*/

	MLOG("Step >> Main Steps Done!\n\n");

	OUTPUT(X, Y, Z, Xst, Xed, Yst, Yed, s, myrank, size, other, x_sec, y_sec, argv[1], local_image, image, rankinfo, nodes);

	arrayFree2DF(image);
	arrayFree2DI(rankinfo);
	arrayFree2DF(local_image);

	arrayFree3DI(flags);
	arrayFree5DF(nodes);
	arrayFree4DI(walls);
	arrayFree3DF(temp_right);
	arrayFree3DF(temp_left);
	arrayFree3DF(temp_down);
	arrayFree3DF(temp_up);
	arrayFree2DF(temp_lu);
	arrayFree2DF(temp_ld);
	arrayFree2DF(temp_ru);
	arrayFree2DF(temp_rd);
	arrayFree3DF(temp_right_send);
	arrayFree3DF(temp_left_send);
	arrayFree3DF(temp_down_send);
	arrayFree3DF(temp_up_send);
	arrayFree2DF(temp_lu_send);
	arrayFree2DF(temp_ld_send);
	arrayFree2DF(temp_ru_send);
	arrayFree2DF(temp_rd_send);

	MLOG("LBM-simulation Done!\n");

  MPI_Finalize();

	return 0;
}
