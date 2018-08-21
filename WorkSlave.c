#include "slave.h"
#include "simd.h"

typedef float Real;

__thread_local const int _dfInv[19] = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};
__thread_local const int _e_x[20]__attribute__((aligned(32))) = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0,0};
__thread_local const int _e_y[20]__attribute__((aligned(32))) = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0,0};
__thread_local const int _e_z[20]__attribute__((aligned(32))) = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0,0};
__thread_local Real _w[20]__attribute__((aligned(32))) = {
 	(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
 	(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
 	(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.),(0.0) };

__thread_local int my_id;

__thread_local Real *****nodes;
__thread_local char ****walls;

__thread_local Real node[4][12][12][19] __attribute__((aligned(128)));
__thread_local Real ans[10][10][19] __attribute__((aligned(128)));
__thread_local char wall[10][10][20] __attribute__((aligned(128)));
__thread_local Real pre[10][10][19] __attribute__((aligned(128)));

__thread_local volatile unsigned long nodeflag[4], ansflag[2], wallflag[2], putflag[2];
__thread_local int cs, os, now;

__thread_local Real rho, u_x, u_y, u_z,nfSub[20],feq[20];
__thread_local floatv4 rxyzV, eCoefV[20];
__thread_local Real Qo, omegaNew, S;
__thread_local Real _nu,_omega,_CSmago;

extern Real omega,CSmago;

inline Real sqr(Real s) { return s * s; }

void get_wall(int j, int i, int k, int l)
{
    wallflag[l] = 0;
    if (l == 0)
    {
        athread_get(PE_MODE, &walls[j * 1250 + i * 5 + k / 10][0][0][0], &wall[0][0][0], 1000, &wallflag[l], 0, 0, 0);
    }
    else
    {
        athread_get(PE_MODE, &walls[j * 1250 + i * 5 + k / 10][5][0][0], &wall[5][0][0], 1000, &wallflag[l], 0, 0, 0);
    }
}
void get_ans(int j, int i, int k, int l)
{
    ansflag[l] = 0;
    if (l == 0)
    {
        athread_get(PE_MODE, &nodes[cs][i][j][k][0], &pre[0][0][0], 3800, &ansflag[0], 0, 4825240, 760);
    }
    else
    {
        athread_get(PE_MODE, &nodes[cs][i + 5][j][k][0], &pre[5][0][0], 3800, &ansflag[1], 0, 4825240, 760);
    }
}
void get_node(int j, int i, int k, int p)
{
    nodeflag[p] = 0;
    if (k == 0)
    {
        athread_get(PE_MODE, &nodes[os][i - 1][j][k][0], &node[p][0][0][0], 10944, &nodeflag[p], 0, 4825088, 912);
    }
    else if (k == 490)
    {
        athread_get(PE_MODE, &nodes[os][i - 1][j][k - 2][0], &node[p][0][0][0], 10944, &nodeflag[p], 0, 4825088, 912);
    }
    else
    {
        athread_get(PE_MODE, &nodes[os][i - 1][j][k - 1][0], &node[p][0][0][0], 10944, &nodeflag[p], 0, 4825088, 912);
    }
}
void collide(Real *nodesPerCell)
{
    int l;
    rxyzV = 0;
    rxyzV += (floatv4)(nodesPerCell[0]) * eCoefV[0];
    rxyzV += (floatv4)(nodesPerCell[1]) * eCoefV[1];
    rxyzV += (floatv4)(nodesPerCell[2]) * eCoefV[2];
    rxyzV += (floatv4)(nodesPerCell[3]) * eCoefV[3];
    rxyzV += (floatv4)(nodesPerCell[4]) * eCoefV[4];
    rxyzV += (floatv4)(nodesPerCell[5]) * eCoefV[5];
    rxyzV += (floatv4)(nodesPerCell[6]) * eCoefV[6];
    rxyzV += (floatv4)(nodesPerCell[7]) * eCoefV[7];
    rxyzV += (floatv4)(nodesPerCell[8]) * eCoefV[8];
    rxyzV += (floatv4)(nodesPerCell[9]) * eCoefV[9];
    rxyzV += (floatv4)(nodesPerCell[10]) * eCoefV[10];
    rxyzV += (floatv4)(nodesPerCell[11]) * eCoefV[11];
    rxyzV += (floatv4)(nodesPerCell[12]) * eCoefV[12];
    rxyzV += (floatv4)(nodesPerCell[13]) * eCoefV[13];
    rxyzV += (floatv4)(nodesPerCell[14]) * eCoefV[14];
    rxyzV += (floatv4)(nodesPerCell[15]) * eCoefV[15];
    rxyzV += (floatv4)(nodesPerCell[16]) * eCoefV[16];
    rxyzV += (floatv4)(nodesPerCell[17]) * eCoefV[17];
    rxyzV += (floatv4)(nodesPerCell[18]) * eCoefV[18];

    rho = simd_vextf0(rxyzV);
    u_x = simd_vextf1(rxyzV) / rho;
    u_y = simd_vextf2(rxyzV) / rho;
    u_z = simd_vextf3(rxyzV) / rho;
    //rho = (u_x * u_x + u_y * u_y + u_z * u_z)
    double tmp2 = 1.5 * (u_x * u_x + u_y * u_y + u_z * u_z);
	for (l = 0; l < 19; l++) {
		const Real tmp = (_e_x[l] * u_x + _e_y[l] * u_y + _e_z[l] * u_z);
		feq[l] = _w[l] * rho * (1.0-tmp2+(3.0*tmp)+(4.5 * tmp * tmp));
        nfSub[l]=nodesPerCell[l]-feq[l];
	}
    Qo=0;
    Qo += sqr(nfSub[2] + nfSub[3] + nfSub[6] + nfSub[7] + nfSub[8] + nfSub[9] + nfSub[14] + nfSub[15] + nfSub[16] + nfSub[17]);
    Qo += sqr(nfSub[6] - nfSub[7] - nfSub[8] + nfSub[9]);
    Qo += sqr(nfSub[14] - nfSub[15] - nfSub[16] + nfSub[17]);
    double sij=nfSub[6]-nfSub[7]-nfSub[8]+nfSub[9];
	Qo += sij * sij;
    Qo += sqr(nfSub[0] + nfSub[1] + nfSub[6] + nfSub[7] + nfSub[8] + nfSub[9] + nfSub[10] + nfSub[11] + nfSub[12] + nfSub[13]);
    Qo += sqr(nfSub[10] - nfSub[11] - nfSub[12] + nfSub[13]);
    Qo += sqr(nfSub[14] - nfSub[15] - nfSub[16] + nfSub[17]);
    Qo += sqr(nfSub[10] - nfSub[11] - nfSub[12] + nfSub[13]);
    Qo += sqr(nfSub[4] + nfSub[5] + nfSub[10] + nfSub[11] + nfSub[12] + nfSub[13] + nfSub[14] + nfSub[15] + nfSub[16] + nfSub[17]);
	S = (-_nu + sqrt(_nu * _nu + 18 * _CSmago * _CSmago * sqrt(Qo))) / 6.0 / _CSmago / _CSmago;
	omegaNew = 1.0 / (3.0 * (_nu + _CSmago * _CSmago * S) + 0.5);
    rho = (1.0-omegaNew);
	nodesPerCell[0] = rho * nodesPerCell[0] +omegaNew * feq[0];
    nodesPerCell[1] = rho * nodesPerCell[1] +omegaNew * feq[1];
    nodesPerCell[2] = rho * nodesPerCell[2] +omegaNew * feq[2];
    nodesPerCell[3] = rho * nodesPerCell[3] +omegaNew * feq[3];
    nodesPerCell[4] = rho * nodesPerCell[4] +omegaNew * feq[4];
    nodesPerCell[5] = rho * nodesPerCell[5] +omegaNew * feq[5];
    nodesPerCell[6] = rho * nodesPerCell[6] +omegaNew * feq[6];
    nodesPerCell[7] = rho * nodesPerCell[7] +omegaNew * feq[7];
    nodesPerCell[8] = rho * nodesPerCell[8] +omegaNew * feq[8];
    nodesPerCell[9] = rho * nodesPerCell[9] +omegaNew * feq[9];
    nodesPerCell[10] = rho * nodesPerCell[10] +omegaNew * feq[10];
    nodesPerCell[11] = rho * nodesPerCell[11] +omegaNew * feq[11];
    nodesPerCell[12] = rho * nodesPerCell[12] +omegaNew * feq[12];
    nodesPerCell[13] = rho * nodesPerCell[13] +omegaNew * feq[13];
    nodesPerCell[14] = rho * nodesPerCell[14] +omegaNew * feq[14];
    nodesPerCell[15] = rho * nodesPerCell[15] +omegaNew * feq[15];
    nodesPerCell[16] = rho * nodesPerCell[16] +omegaNew * feq[16];
    nodesPerCell[17] = rho * nodesPerCell[17] +omegaNew * feq[17];
    nodesPerCell[18] = rho * nodesPerCell[18] +omegaNew * feq[18];
}
void workParallel(long *t)
{
    my_id = athread_get_id(-1);
    if (my_id >= 60)
        return;

    nodes = t[0];
    walls = t[1];
    cs = t[2];
    os = cs ^ 1;
    _omega = omega;
	_CSmago = CSmago;
	_nu = (2.0 / _omega - 1.0) / 6.0;

    int i, j, k, l, p, q;
    int inv, tq;

    int z = my_id % 5;
    int y = (my_id / 5) % 4;
    int x = my_id / 20;
    int ze = 100;
    int ye = (y == 3) ? 32 : 31;
    int xe = (x == 2) ? 90 : 80;
    z *= 100;
    y *= 31;
    x *= 80;
    putflag[0] = 1;
    putflag[1] = 1;

    for (i = 0; i < 20; ++i)
        eCoefV[i] = simd_set_floatv4(1.0, _e_x[i], _e_y[i], _e_z[i]);

    for (i = 1; i <= xe; i += 10)
    {
        for (k = 0; k < ze; k += 10)
        {
            now = 0;
            get_node(y, i + x, z + k, 0);
            get_node(y + 1, i + x, z + k, 1);
            get_node(y + 2, i + x, z + k, 2);
            get_wall(y, i + x - 1, z + k, 0);
            get_ans(y + 1, i + x, z + k, 0);
            while (nodeflag[0] != 1);
            while (nodeflag[1] != 1);
            if (z + k == 0)
                tq = 0;
            else if (z + k == 490)
                tq = 2;
            else
                tq = 1;
            for (j = 1; j <= ye; j++)
            {
                now++;
                if (now > 3)
                    now -= 4;
                if (j < ye)
                    get_node(j + y + 2, i + x, z + k, now ^ 2);
                get_wall(j + y - 1, i + x - 1, z + k, 1);
                get_ans(j + y, i + x, z + k, 1);
                while (wallflag[0] != 1);
                while (ansflag[0] != 1);
                while (putflag[0] != 1);
                while (nodeflag[(now + 1) & 3] != 1);
                for (p = 0; p < 5; p++)
                    for (q = 0; q < 10; q++)
                    {
                        if (wall[p][q][19] == 0)
                        {
                            ans[p][q][ 0] = node[(now + 3) & 3][p + 1][q + tq    ][ 0];
                            ans[p][q][ 1] = node[(now + 5) & 3][p + 1][q + tq    ][ 1];
                            ans[p][q][ 2] = node[(now + 4) & 3][p    ][q + tq    ][ 2];
                            ans[p][q][ 3] = node[(now + 4) & 3][p + 2][q + tq    ][ 3];
                            ans[p][q][ 4] = node[(now + 4) & 3][p + 1][q + tq - 1][ 4];
                            ans[p][q][ 5] = node[(now + 4) & 3][p + 1][q + tq + 1][ 5];
                            ans[p][q][ 6] = node[(now + 3) & 3][p    ][q + tq    ][ 6];
                            ans[p][q][ 7] = node[(now + 3) & 3][p + 2][q + tq    ][ 7];
                            ans[p][q][ 8] = node[(now + 5) & 3][p    ][q + tq    ][ 8];
                            ans[p][q][ 9] = node[(now + 5) & 3][p + 2][q + tq    ][ 9];
                            ans[p][q][10] = node[(now + 3) & 3][p + 1][q + tq - 1][10];
                            ans[p][q][11] = node[(now + 3) & 3][p + 1][q + tq + 1][11];
                            ans[p][q][12] = node[(now + 5) & 3][p + 1][q + tq - 1][12];
                            ans[p][q][13] = node[(now + 5) & 3][p + 1][q + tq + 1][13];
                            ans[p][q][14] = node[(now + 4) & 3][p    ][q + tq - 1][14];
                            ans[p][q][15] = node[(now + 4) & 3][p    ][q + tq + 1][15];
                            ans[p][q][16] = node[(now + 4) & 3][p + 2][q + tq - 1][16];
                            ans[p][q][17] = node[(now + 4) & 3][p + 2][q + tq + 1][17];
                            ans[p][q][18] = node[(now + 4) & 3][p + 1][q + tq    ][18];
                            collide(ans[p][q]);
                        }
                        else if (wall[p][q][19] == 3)
                        {
                            for (l = 0; l < 19; l++)
                            {
                                inv = _dfInv[l];
                                if (wall[p][q][l])
                                {
                                    ans[p][q][l] = node[now][p + 1][q + tq][inv];
                                }
                                else
                                {
                                    ans[p][q][l] = node[(now + _e_y[inv] + 4) & 3][p + 1 + _e_x[inv]][q + tq + _e_z[inv]][l];
                                }
                            }
                            collide(ans[p][q]);
                        }
                        else
                            memcpy(&ans[p][q][0], &pre[p][q][0], 76);
                    }
                putflag[0] = 0;
                athread_put(PE_MODE, &ans[0][0][0], &nodes[cs][i + x][j + y][z + k][0], 3800, &putflag[0], 4825240, 760);
                if (j < ye)
                {
                    get_wall(j + y, i + x - 1, z + k, 0);
                    get_ans(j + y + 1, i + x, z + k, 0);
                }
                while (wallflag[1] != 1);
                while (ansflag[1] != 1);

                for (p = 5; p < 10; p++)
                    for (q = 0; q < 10; q++)
                    {
                        if (wall[p][q][19] == 0)
                        {
                            ans[p][q][ 0] = node[(now + 3) & 3][p + 1][q + tq    ][ 0];
                            ans[p][q][ 1] = node[(now + 5) & 3][p + 1][q + tq    ][ 1];
                            ans[p][q][ 2] = node[(now + 4) & 3][p    ][q + tq    ][ 2];
                            ans[p][q][ 3] = node[(now + 4) & 3][p + 2][q + tq    ][ 3];
                            ans[p][q][ 4] = node[(now + 4) & 3][p + 1][q + tq - 1][ 4];
                            ans[p][q][ 5] = node[(now + 4) & 3][p + 1][q + tq + 1][ 5];
                            ans[p][q][ 6] = node[(now + 3) & 3][p    ][q + tq    ][ 6];
                            ans[p][q][ 7] = node[(now + 3) & 3][p + 2][q + tq    ][ 7];
                            ans[p][q][ 8] = node[(now + 5) & 3][p    ][q + tq    ][ 8];
                            ans[p][q][ 9] = node[(now + 5) & 3][p + 2][q + tq    ][ 9];
                            ans[p][q][10] = node[(now + 3) & 3][p + 1][q + tq - 1][10];
                            ans[p][q][11] = node[(now + 3) & 3][p + 1][q + tq + 1][11];
                            ans[p][q][12] = node[(now + 5) & 3][p + 1][q + tq - 1][12];
                            ans[p][q][13] = node[(now + 5) & 3][p + 1][q + tq + 1][13];
                            ans[p][q][14] = node[(now + 4) & 3][p    ][q + tq - 1][14];
                            ans[p][q][15] = node[(now + 4) & 3][p    ][q + tq + 1][15];
                            ans[p][q][16] = node[(now + 4) & 3][p + 2][q + tq - 1][16];
                            ans[p][q][17] = node[(now + 4) & 3][p + 2][q + tq + 1][17];
                            ans[p][q][18] = node[(now + 4) & 3][p + 1][q + tq    ][18];
                            collide(ans[p][q]);
                        }
                        else if (wall[p][q][19] == 3)
                        {
                            for (l = 0; l < 19; l++)
                            {
                                inv = _dfInv[l];
                                if (wall[p][q][l])
                                {
                                    ans[p][q][l] = node[now][p + 1][q + tq][inv];
                                }
                                else
                                {
                                    ans[p][q][l] = node[(now + _e_y[inv] + 4) & 3][p + 1 + _e_x[inv]][q + tq + _e_z[inv]][l];
                                }
                            }
                            collide(ans[p][q]);
                        }
                        else
                            memcpy(&ans[p][q][0], &pre[p][q][0], 76);
                    }
                putflag[1] = 0;
                athread_put(PE_MODE, &ans[5][0][0], &nodes[cs][i + x + 5][j + y][z + k][0], 3800, &putflag[1], 4825240, 760);
            }
        }
    }
}
