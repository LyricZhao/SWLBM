# include "Argument.h"
# include "myArgument.h"

# include "simd.h"
# include <mpi.h>
# include <math.h>

Real rho, u_x, u_y, u_z, uxyzConst, nf6, nf10, nf14, nfSub[20];

floatv4 tmpV, fiV, rxyzV, eCoefV[20], feqV[5];
floatv4 *e_xV = (floatv4*) e_xM, *e_yV = (floatv4*) e_yM, *e_zV = (floatv4*) e_zM, *wV = (floatv4*) wM;
floatv4 const1V = 1.0, const3V = 3.0, const4p5V = 4.5;
floatv4 uxyzConstV, u_xV, u_yV, u_zV, omegaNewV, rhoConstV;
floatv4 *nfSubV = (floatv4*) nfSub;

Real Qo, omegaNew, S, *nodesPerCell, *feq = (Real*) feqV;
Real *npc0, ****npc1;

inline Real sqr(Real s) { return s * s; }

void computeOneStep(Real *****nodes, int ****walls, int ***flags, int Xst, int Xed, int Yst, int Yed, int nz, int current) {
	int i, j, k, l, m, inv, other = current ^ 1;

	// Prepare
	for(i = 0; i < 20; ++ i) {
		eCoefV[i] = simd_set_floatv4(1.0, e_xM[i], e_yM[i], e_zM[i]);
	}

	// Main Computing Part
	for(i = 1; i < Xed - Xst + 1; ++ i) {
		for(j = 1; j < Yed - Yst + 1; ++ j) {

			int* flagsPerRow = flags[i][j];

			for(k = 0; k < nz; ++ k) {
        if(flagsPerRow[k] == FLUID) {
          npc1 = nodes[other];
          npc0 = nodes[current][i][j][k];

          npc0[ 0] = npc1[i    ][j - 1][k    ][ 0];
          npc0[ 1] = npc1[i    ][j + 1][k    ][ 1];
          npc0[ 2] = npc1[i - 1][j    ][k    ][ 2];
          npc0[ 3] = npc1[i + 1][j    ][k    ][ 3];
          npc0[ 4] = npc1[i    ][j    ][k - 1][ 4];
          npc0[ 5] = npc1[i    ][j    ][k + 1][ 5];
          npc0[ 6] = npc1[i - 1][j - 1][k    ][ 6];
          npc0[ 7] = npc1[i + 1][j - 1][k    ][ 7];
          npc0[ 8] = npc1[i - 1][j + 1][k    ][ 8];
          npc0[ 9] = npc1[i + 1][j + 1][k    ][ 9];
          npc0[10] = npc1[i    ][j - 1][k - 1][10];
          npc0[11] = npc1[i    ][j - 1][k + 1][11];
          npc0[12] = npc1[i    ][j + 1][k - 1][12];
          npc0[13] = npc1[i    ][j + 1][k + 1][13];
          npc0[14] = npc1[i - 1][j    ][k - 1][14];
          npc0[15] = npc1[i - 1][j    ][k + 1][15];
          npc0[16] = npc1[i + 1][j    ][k - 1][16];
          npc0[17] = npc1[i + 1][j    ][k + 1][17];
          npc0[18] = npc1[i    ][j    ][k    ][18];

        } else if(flags[i][j][k] == BOUNCE) {
          for(l = 0; l < 19; ++ l) {
            inv = dfInv[l];
            if(walls[i - 1][j - 1][k][l]) {
              nodes[current][i][j][k][l] = nodes[other][i][j][k][inv];
            } else {
              nodes[current][i][j][k][l] = nodes[other][i + e_x[inv]][j + e_y[inv]][k + e_z[inv]][l];
            }
          }
        }

				if(flagsPerRow[k] == FLUID || flagsPerRow[k] == BOUNCE) {

					rxyzV = 0;
					nodesPerCell = nodes[current][i][j][k];
					for(l = 0; l < 19; ++ l) {
						rxyzV += (floatv4)(nodesPerCell[l]) * eCoefV[l];
					} // Loop l, Vectorized
					rho = simd_vextf0(rxyzV);
					u_x = simd_vextf1(rxyzV) / rho;
					u_y = simd_vextf2(rxyzV) / rho;
					u_z = simd_vextf3(rxyzV) / rho;

					uxyzConst = 1.5 * (u_x * u_x + u_y * u_y + u_z * u_z);
					uxyzConstV = (floatv4) uxyzConst, rhoConstV = (floatv4) rho;
					u_xV = u_x, u_yV = u_y, u_zV = u_z;
					for(m = 0; m < 5; ++ m) {
						tmpV = e_xV[m] * u_xV + e_yV[m] * u_yV + e_zV[m] * u_zV;
						feqV[m] = wV[m] * rhoConstV * (const1V - uxyzConstV + const3V * tmpV + const4p5V * tmpV * tmpV);
					}  // Loop m, Vectorized

					Qo = 0;
					for(m = 0; m < 4; ++ m) {
						simd_loadu(tmpV, nodesPerCell + (m << 2));
						nfSubV[m] = tmpV - feqV[m];
					} // Loop m, Vectorized
					nfSub[16] = nodesPerCell[16] - feq[16];
					nfSub[17] = nodesPerCell[17] - feq[17];
					nfSub[18] = nodesPerCell[18] - feq[18];

					nf6 = nfSub[6] + nfSub[7] + nfSub[8] + nfSub[9];
					nf10 = nfSub[10] + nfSub[11] + nfSub[12] + nfSub[13];
					nf14 = nfSub[14] + nfSub[15] + nfSub[16] + nfSub[17];
					/* xx */ Qo += sqr(nfSub[2] + nfSub[3] + nf6 + nf14);
					/* yy */ Qo += sqr(nfSub[0] + nfSub[1] + nf6 + nf10);
					/* zz */ Qo += sqr(nfSub[4] + nfSub[5] + nf10 + nf14);
					/* xy */ Qo += 2 * sqr(nfSub[6] - nfSub[7] - nfSub[8] + nfSub[9]);
					/* xz */ Qo += 2 * sqr(nfSub[14] - nfSub[15] - nfSub[16] + nfSub[17]);
					/* yz */ Qo += 2 * sqr(nfSub[10] - nfSub[11] - nfSub[12] + nfSub[13]);

					nu = (2.0 / omega - 1.0) / 6.0;
					S = (-nu + sqrt(nu * nu + 18 * CSmago * CSmago * sqrt(Qo))) / 6.0;
					omegaNew = 1.0 / (3.0 * (nu + S) + 0.5);

					omegaNewV = (floatv4) omegaNew;
					for(m = 0; m < 4; ++ m) {
						simd_loadu(tmpV, nodesPerCell + (m << 2));
						tmpV = (const1V - omegaNewV) * tmpV + omegaNewV * feqV[m];
						simd_storeu(tmpV, nodesPerCell + (m << 2));
					} // Loop m, Vectorized
					nodesPerCell[16] = (1.0 - omegaNew) * nodesPerCell[16] + omegaNew * feq[16];
					nodesPerCell[17] = (1.0 - omegaNew) * nodesPerCell[17] + omegaNew * feq[17];
					nodesPerCell[18] = (1.0 - omegaNew) * nodesPerCell[18] + omegaNew * feq[18];

				} // If flag
			} // Loop k
		} // Loop j
	} // Loop i
}
