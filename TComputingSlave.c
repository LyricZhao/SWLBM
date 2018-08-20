# include "simd.h"
# include "slave.h"

# include <math.h>
# include <assert.h>

# include "Argument.h"
# include "ArgumentMT.h"

/* Parameters */
__thread_kernel("ldm") Real *****nodes;
__thread_kernel("ldm") int ****walls, ***flags, Xst, Xed, Yst, Yed, nz, current, other;

/* Thread */
__thread_kernel("ldm") volatile long wait = 0;
__thread_kernel("ldm") int threadID, threadST, threadED, statusI[3];
__thread_kernel("ldm") Real nodesCurrentLocal[MAXJ - 2][MAXK][MAXL] __attribute__((aligned(128)));

/* Main Vars */
__thread_kernel("ldm") Real nodesOtherLocal[MAXI][MAXJ][MAXK + 2][MAXL] __attribute__((aligned(128)));
__thread_kernel("ldm") int wallsLocal[MAXJ - 2][MAXK][MAXL] __attribute__((aligned(128)));
__thread_kernel("ldm") int flagsLocal[MAXJ - 2][MAXK] __attribute__((aligned(128)));

/* Collide */
__thread_kernel("ldm") Real Qo, S, omegaNew, rho, u_x, u_y, u_z, uxyzConst, nf6, nf10, nf14, *feq;
__thread_kernel("ldm") Real nfSub[20] __attribute__((aligned(16)));

/* Vector */
__thread_kernel("ldm") floatv4 tmpV, fiV, rxyzV;
__thread_kernel("ldm") floatv4 uxyzConstV, u_xV, u_yV, u_zV, omegaNewV, rhoConstV;
__thread_kernel("ldm") floatv4 const1V = 1.0, const3V = 3.0, const4p5V = 4.5;
__thread_kernel("ldm") floatv4 *e_xV = (floatv4*) e_xM, *e_yV = (floatv4*) e_yM, *e_zV = (floatv4*) e_zM, *wV = (floatv4*) wM, *nfSubV;
__thread_kernel("ldm") floatv4 eCoefV[20] __attribute__((aligned(16)));
__thread_kernel("ldm") floatv4 feqV[5] __attribute__((aligned(16)));

inline void setParas(long *para) {
  nodes = (Real *****) para[0]; walls = (int ****) para[1]; flags = (int ***) para[2];
  Xst = para[3]; Xed = para[4]; Yst = para[5];
  Yed = para[6]; nz = para[7]; current = para[8];
  other = current ^ 1, nfSubV = (floatv4 *) nfSub, feq = (Real *) feqV;
}

inline void threadInit() {
  threadID = athread_get_id(-1);
  int length = Yed - Yst;
  threadST = BLOCK_LOW(threadID, MT_NUMS, length) + 1;
  threadED = BLOCK_HIH(threadID, MT_NUMS, length) + 2;
  /* NOTE: [threadST, threadED) */
  assert((threadED - threadST) <= (MAXJ - 2));
  statusI[0] = statusI[1] = statusI[2] = 0;
}

inline Real sqr(Real s) {
  return s * s;
}

inline void DMA_Get(void *source, void *target, int size) {
  wait = 0;
  athread_get(PE_MODE, source, target, size, (void *) &wait, 0, 0, 0);
  while(!wait);
}

inline void DMA_Put(void *source, void *target, int size) {
  wait = 0;
  athread_put(PE_MODE, source, target, size, (void *) &wait, 0, 0);
  while(!wait);
}

/* Input: 250 * 125 * 500 */
void computeOneStepParallel(long *para) {

  setParas(para);
  threadInit();

  int i, j, k, l, t, m, kST, inv;
  for(kST = 0; kST < nz; kST += MAXK) {
    for(i = 1; i < Xed - Xst + 1; ++ i) {

      /* Communicate */
      for(t = i - 1; t < i + 2; ++ t) {
        if(statusI[t % 3] == t) continue;
        statusI[t % 3] = t;
        for(j = threadST - 1; j < threadED + 1; ++ j) {
          if(kST == 0) {
            DMA_Get(&nodes[other][t][j][kST][0], &nodesOtherLocal[t % 3][j & 3][1][0], (MAXK + 1) * MAXL * sizeof(Real));
          }
          else if(kST + MAXK == nz) {
            DMA_Get(&nodes[other][t][j][kST][0], &nodesOtherLocal[t % 3][j & 3][0][0], (MAXK + 1) * MAXL * sizeof(Real));
          }
          else {
            DMA_Get(&nodes[other][t][j][kST][0], &nodesOtherLocal[t % 3][j & 3][0][0], (MAXK + 2) * MAXL * sizeof(Real));
          }
        }
      }
      for(j = threadST; j < threadED; ++ j) {
        DMA_Get(&nodes[current][i][j][kST][0], &nodesCurrentLocal[j & 1][0][0], MAXK * MAXL * sizeof(Real));
        DMA_Get(&walls[i - 1][j - 1][kST][0], &wallsLocal[j & 1][0][0], MAXK * MAXL * sizeof(int));
        DMA_Get(&flags[i][j][kST], &flagsLocal[j & 1][0], MAXK * sizeof(int));
      }

      /* Compute */
      for(j = threadST; j < threadED; ++ j) {
        int *flagsM = (int *) flagsLocal[j & 1], **wallsM = (int **) wallsLocal[j & 1];
        for(k = 0; k < MAXK; ++ k) {
          Real *npc0 = nodesCurrentLocal[j & 1][k];

          /* Stream */
          if(flagsM[k] == FLUID) {
            int is1 = (i - 1) % 3, in0 = i % 3, ia1 = (i + 1) % 3;
            int js1 = (j - 1) & 3, jn0 = j & 3, ja1 = (j + 1) & 3;
            npc0[ 0] = nodesOtherLocal[in0][js1][k + 1][ 0];
            npc0[ 1] = nodesOtherLocal[in0][ja1][k + 1][ 1];
            npc0[ 2] = nodesOtherLocal[is1][jn0][k + 1][ 2];
            npc0[ 3] = nodesOtherLocal[ia1][jn0][k + 1][ 3];
            npc0[ 4] = nodesOtherLocal[in0][jn0][k    ][ 4];
            npc0[ 5] = nodesOtherLocal[in0][jn0][k + 2][ 5];
            npc0[ 6] = nodesOtherLocal[is1][js1][k + 1][ 6];
            npc0[ 7] = nodesOtherLocal[ia1][js1][k + 1][ 7];
            npc0[ 8] = nodesOtherLocal[is1][ja1][k + 1][ 8];
            npc0[ 9] = nodesOtherLocal[ia1][ja1][k + 1][ 9];
            npc0[10] = nodesOtherLocal[in0][js1][k    ][10];
            npc0[11] = nodesOtherLocal[in0][js1][k + 2][11];
            npc0[12] = nodesOtherLocal[in0][ja1][k    ][12];
            npc0[13] = nodesOtherLocal[in0][ja1][k + 2][13];
            npc0[14] = nodesOtherLocal[is1][jn0][k    ][14];
            npc0[15] = nodesOtherLocal[is1][jn0][k + 2][15];
            npc0[16] = nodesOtherLocal[ia1][jn0][k    ][16];
            npc0[17] = nodesOtherLocal[ia1][jn0][k + 2][17];
            npc0[18] = nodesOtherLocal[in0][jn0][k + 1][18];
          } else if(flagsM[k] == BOUNCE) {
            for(l = 0; l < 19; ++ l) {
              inv = dfInv[l];
              if(wallsM[k][l]) {
                npc0[l] = nodesOtherLocal[i % 3][j & 3][k + 1][inv];
              } else {
                npc0[l] = nodesOtherLocal[(i + e_x[inv]) % 3][(j + e_y[inv]) & 3][k + e_z[inv] + 1][l];
              }
            }
          }

          /* Collide */
          if(flagsM[k] == FLUID || flagsM[k] == BOUNCE) {
  					rxyzV = 0;
  					for(l = 0; l < 19; ++ l) {
  						rxyzV += (floatv4)(npc0[l]) * eCoefV[l];
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
  						simd_loadu(tmpV, npc0 + (m << 2));
  						nfSubV[m] = tmpV - feqV[m];
  					} // Loop m, Vectorized
  					nfSub[16] = npc0[16] - feq[16];
  					nfSub[17] = npc0[17] - feq[17];
  					nfSub[18] = npc0[18] - feq[18];

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
  						simd_loadu(tmpV, npc0 + (m << 2));
  						tmpV = (const1V - omegaNewV) * tmpV + omegaNewV * feqV[m];
  						simd_storeu(tmpV, npc0 + (m << 2));
  					} // Loop m, Vectorized
  					npc0[16] = (1.0 - omegaNew) * npc0[16] + omegaNew * feq[16];
  					npc0[17] = (1.0 - omegaNew) * npc0[17] + omegaNew * feq[17];
  					npc0[18] = (1.0 - omegaNew) * npc0[18] + omegaNew * feq[18];
  				}
        } // Loop k
        DMA_Put(&nodesCurrentLocal[j & 1][0][0], &nodes[current][i][j][kST][0], MAXK * MAXL * sizeof(Real));
      }
    }
  }
}
