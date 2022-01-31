//ガウス変数の生成
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include  "MT.h"
#define N 1e5 //サンプル数

void fwht(int, int, int, double, double*, double*, double*, double*);
/* n = (int) pow(2.0, (double)m);
   nh = (double)n / 2.0;
   rrn = 1.0 / sqrt(double(n));
   dim xI[n], xR[n], yR[n], yI[n]; */

void fwht(int m, int n, int nh, double rrn, double* xR, double* xI, double* yR, double* yI) {

    //変数宣言
    int i, i0, i1, j;

    double* zR, * zI;
    zR = calloc(N, sizeof(double));
    if (zR == NULL) {
        puts("配列zR[n]の記憶領域の確保に失敗");
        exit(1);
    }

    zI = calloc(N, sizeof(double));
    if (zI == NULL) {
        puts("配列zI[n]の記憶領域の確保に失敗");
        exit(1);
    }

    //入力データ格納
    for (i = 0; i < n; i++) {
        zR[i] = xR[i];  zI[i] = xI[i];
    }

    //高速ウォルシュアダマール変換(FWHT)
    for (j = 0; j < m; j++) {
        i0 = 0;

        for (i = 0; (double)i < nh; i++) {
            i1 = i0 + 1;
            yR[i] = zR[i0] + zR[i1];  yI[i] = zI[i0] + zI[i1];
            yR[i + nh] = zR[i0] - zR[i1];  yI[i + nh] = zI[i0] - zI[i1];
            i0 = i1 + 1;
        }

        for (i = 0; i < n; i++) {
            zR[i] = yR[i];  zI[i] = yI[i];
        }

    }

    for (i = 0; i < n; i++) {
        yR[i] = rrn * zR[i];  yI[i] = rrn * zI[i];
    }

    free(zR); free(zI);

    return;
}
