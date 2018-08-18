# include "Argument.h"
# include "myArgument.h"

# include "simd.h"

Real *npc0, ****npc1;

void stream(Real *****nodes, int ****walls, int ***flags, int Xst, int Xed, int Yst, int Yed, int nz, int current, int other) {
  int i, j, k, l, inv;

  for(i = 1; i < Xed - Xst + 1; ++ i) {
		for(j = 1; j < Yed - Yst + 1; ++ j) {
			for(k = 0; k < nz; ++ k) {
				if(flags[i][j][k] == FLUID) {
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
			}
		}
	}
}
