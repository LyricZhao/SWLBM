# include "simd.h"
# include "slave.h"

# include <math.h>
# include <stdio.h>
# include <assert.h>

# include "ArgumentMT.h"

/* Parameters */
__thread_local_fix Real *****nodes;
__thread_local_fix char ****walls;
__thread_local_fix int Xst, Xed, Yst, Yed, nz, current, other;

/* Thread */
__thread_local_fix volatile long wait, target, putFlag, putTarget, waitWalls, targetWalls;
__thread_local_fix int threadID, threadST, threadED, statusI[3];
__thread_local_fix Real nodesCurrentLocal[MAXJ - 2][MAXK][MAXL] __attribute__((aligned(128)));

/* Main Vars */
__thread_local_fix Real nodesOtherLocal[MAXI][MAXJ][MAXK + 2][MAXL] __attribute__((aligned(128)));
__thread_local_fix char wallsLocal[MAXJ - 2][MAXK][MAXL + 1] __attribute__((aligned(128)));

/* Collide */
__thread_local_fix Real _nu, _omega, _CSmago;
__thread_local_fix Real Qo, S, omegaNew, rho, u_x, u_y, u_z, uxyzConst, nf6, nf10, nf14, *npc0;
__thread_local_fix Real nfSub[20] __attribute__((aligned(32)));
__thread_local_fix Real feq[20] __attribute__((aligned(32)));

/* Vector */
__thread_local_fix floatv4 tmpV, fiV, rxyzV;
__thread_local_fix floatv4 uxyzConstV, u_xV, u_yV, u_zV, omegaNewV, rhoConstV;
__thread_local_fix floatv4 const1V = 1.0, const3V = 3.0, const4p5V = 4.5;
__thread_local_fix floatv4 eCoefV[20] __attribute__((aligned(32)));

/* Extern */
extern Real nu, omega, CSmago;

inline void setParas(long *para) {
  nodes = (Real *****) para[0]; walls = (char ****) para[1];
  Xst = para[3]; Xed = para[4]; Yst = para[5];
  Yed = para[6]; nz = para[7]; current = para[8];
  other = current ^ 1;
}

inline void threadInit() {
  threadID = athread_get_id(-1);
  int length = Yed - Yst;
  threadST = BLOCK_LOW(threadID, MT_NUMS, length) + 1;
  threadED = BLOCK_HIH(threadID, MT_NUMS, length) + 2;
  /* NOTE: [threadST, threadED) */
  assert((threadED - threadST) <= (MAXJ - 2));
  statusI[0] = statusI[1] = statusI[2] = -1;
}

inline Real sqr(Real s) {
  return s * s;
}

void CommunicateGetWithoutWalls(int i, int kST) {
  int j, t, p = i & 1;
  wait = target = 0;
  for(t = i - 1; t < i + 2; ++ t) {
    if(statusI[t % 3] == t) continue;
    statusI[t % 3] = t;
    for(j = threadST - 1; j < threadED + 1; ++ j) {
      if(kST == 0) {
        ++ target;
        athread_get(PE_MODE, &nodes[other][t][j][kST][0], &nodesOtherLocal[t % 3][j & 3][1][0], (MAXK + 1) * MAXL * sizeof(Real), (void*) &wait, 0, 0, 0);
      }
      else if(kST + MAXK == nz) {
        ++ target;
        athread_get(PE_MODE, &nodes[other][t][j][kST - 1][0], &nodesOtherLocal[t % 3][j & 3][0][0], (MAXK + 1) * MAXL * sizeof(Real), (void*) &wait, 0, 0, 0);
      }
      else {
        ++ target;
        athread_get(PE_MODE, &nodes[other][t][j][kST - 1][0], &nodesOtherLocal[t % 3][j & 3][0][0], (MAXK + 2) * MAXL * sizeof(Real), (void*) &wait, 0, 0, 0);
      }
    }
  }
}

void collideIn() {
  int l;

  rxyzV = 0;
  for (l = 0; l < 19; ++l) {
    rxyzV += (floatv4)(npc0[l]) * eCoefV[l];
  }

  rho = simd_vextf0(rxyzV);
  u_x = simd_vextf1(rxyzV) / rho;
  u_y = simd_vextf2(rxyzV) / rho;
  u_z = simd_vextf3(rxyzV) / rho;

  for (l = 0; l < 19; l++) {
		const Real tmp = (e_xI[l] * u_x + e_yI[l] * u_y + e_zI[l] * u_z);
		feq[l] = wM[l] * rho * (1.0 - (1.5 * (u_x * u_x + u_y * u_y + u_z * u_z)) + (3.0 * tmp) + (4.5 * tmp * tmp));
    nfSub[l] = npc0[l] - feq[l];
	}

  Qo = 0;
  Qo += sqr(nfSub[2] + nfSub[3] + nfSub[6] + nfSub[7] + nfSub[8] + nfSub[9] + nfSub[14] + nfSub[15] + nfSub[16] + nfSub[17]);
  Qo += sqr(nfSub[6] - nfSub[7] - nfSub[8] + nfSub[9]);
  Qo += sqr(nfSub[14] - nfSub[15] - nfSub[16] + nfSub[17]);

  Real sij = nfSub[6] - nfSub[7] - nfSub[8];
  sij += e_yI[9] * e_xI[9] * nfSub[9];
  Qo += sij * sij;

  Qo += sqr(nfSub[0] + nfSub[1] + nfSub[6] + nfSub[7] + nfSub[8] + nfSub[9] + nfSub[10] + nfSub[11] + nfSub[12] + nfSub[13]);
  Qo += sqr(nfSub[10] - nfSub[11] - nfSub[12] + nfSub[13]);
  Qo += sqr(nfSub[14] - nfSub[15] - nfSub[16] + nfSub[17]);
  Qo += sqr(nfSub[10] - nfSub[11] - nfSub[12] + nfSub[13]);
  Qo += sqr(nfSub[4] + nfSub[5] + nfSub[10] + nfSub[11] + nfSub[12] + nfSub[13] + nfSub[14] + nfSub[15] + nfSub[16] + nfSub[17]);
  S = (-_nu + sqrt(_nu * _nu + 18 * _CSmago * _CSmago * sqrt(Qo))) / 6.0 / _CSmago / _CSmago;
  omegaNew = 1.0 / (3.0 * (_nu + _CSmago * _CSmago * S) + 0.5);

	for (l = 0; l < 19; l++)
		npc0[l] = (1.0 - omegaNew) * npc0[l] + omegaNew * feq[l];
}

/* Input: 250 * 125 * 500 */
void computeOneStepParallel(long *para) {

  setParas(para);
  threadInit();

  int i, j, k, l, t, m, kST, inv;
  for(i = 0; i < 20; ++ i) {
    eCoefV[i] = simd_set_floatv4(1.0, e_xM[i], e_yM[i], e_zM[i]);
  }

  _omega = omega;
	_CSmago = CSmago;
	_nu = (2.0 / _omega - 1.0) / 6.0;

  putTarget = putFlag = 0;

  for(kST = 0; kST < nz; kST += MAXK) {

    CommunicateGetWithoutWalls(1, kST);

    for(i = 1; i < Xed - Xst + 1; ++ i) {

      targetWalls = waitWalls = 0;
      for(j = threadST; j < threadED; ++ j) {
        ++ targetWalls;
        athread_get(PE_MODE, &walls[i][j][kST][0], &wallsLocal[j & 1][0][0], MAXK * (MAXL + 1) * sizeof(char), (void*) &waitWalls, 0, 0, 0);
      }
      while(target != wait);
      while(targetWalls != waitWalls);

      /* Compute */
      for(j = threadST; j < threadED; ++ j) {
        for(k = 0; k < MAXK; ++ k) {
          npc0 = nodesCurrentLocal[j & 1][k];

          /* Stream */
          int is1 = (i - 1) % 3, in0 = i % 3, ia1 = (i + 1) % 3;
          int js1 = (j - 1) & 3, jn0 = j & 3, ja1 = (j + 1) & 3;
          if(wallsLocal[j & 1][k][19] == FLUID) {
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
          } else if(wallsLocal[j & 1][k][19] == BOUNCE) {
            for(l = 0; l < 19; ++ l) {
              inv = dfInvM[l];
              if(wallsLocal[j & 1][k][l]) {
                npc0[l] = nodesOtherLocal[i % 3][j & 3][k + 1][inv];
              } else {
                npc0[l] = nodesOtherLocal[(i + e_xI[inv]) % 3][(j + e_yI[inv]) & 3][k + e_zI[inv] + 1][l];
              }
            }
          } else {
            npc0[ 0] = nodesOtherLocal[in0][jn0][k + 1][ 0];
            npc0[ 1] = nodesOtherLocal[in0][jn0][k + 1][ 1];
            npc0[ 2] = nodesOtherLocal[in0][jn0][k + 1][ 2];
            npc0[ 3] = nodesOtherLocal[in0][jn0][k + 1][ 3];
            npc0[ 4] = nodesOtherLocal[in0][jn0][k + 1][ 4];
            npc0[ 5] = nodesOtherLocal[in0][jn0][k + 1][ 5];
            npc0[ 6] = nodesOtherLocal[in0][jn0][k + 1][ 6];
            npc0[ 7] = nodesOtherLocal[in0][jn0][k + 1][ 7];
            npc0[ 8] = nodesOtherLocal[in0][jn0][k + 1][ 8];
            npc0[ 9] = nodesOtherLocal[in0][jn0][k + 1][ 9];
            npc0[10] = nodesOtherLocal[in0][jn0][k + 1][10];
            npc0[11] = nodesOtherLocal[in0][jn0][k + 1][11];
            npc0[12] = nodesOtherLocal[in0][jn0][k + 1][12];
            npc0[13] = nodesOtherLocal[in0][jn0][k + 1][13];
            npc0[14] = nodesOtherLocal[in0][jn0][k + 1][14];
            npc0[15] = nodesOtherLocal[in0][jn0][k + 1][15];
            npc0[16] = nodesOtherLocal[in0][jn0][k + 1][16];
            npc0[17] = nodesOtherLocal[in0][jn0][k + 1][17];
            npc0[18] = nodesOtherLocal[in0][jn0][k + 1][18];
          }
        } // Loop k
      }

      if(i != Xed - Xst) CommunicateGetWithoutWalls(i + 1, kST);

      for(j = threadST; j < threadED; ++ j) {
        for(k = 0; k < MAXK; ++ k) {
          if(wallsLocal[j & 1][k][19] == FLUID || wallsLocal[j & 1][k][19] == BOUNCE) {
              npc0 = nodesCurrentLocal[j & 1][k];
              collideIn();
          }
        }
        ++ putTarget;
        athread_put(PE_MODE, &nodesCurrentLocal[j & 1][0][0], &nodes[current][i][j][kST][0], MAXK * MAXL * sizeof(Real), (void *) &putFlag, 0, 0);
      }
    }
  }
  while(putTarget != putFlag);
}
