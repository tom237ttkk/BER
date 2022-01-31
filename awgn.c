//ガウス変数の生成
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include  "MT.h"

#define GENRAND_MAX 0xffffffffL //MT法により生成される乱数の最大値4294967295
#define gnrnd() ((1.0/(GENRAND_MAX+1.0))*genrand_int32()) // 乱数生成マクロ

void awgn(double, double*, double*);

void awgn(double sigma, double* In, double* Qu) {

    //変数宣言
    double u1, u2, r, theta; //乱数生成用変数

    //(手順1) 区間[0,1)上の一様乱数の生成
    u1 = gnrnd();
    u2 = gnrnd();

    /*(手順2) r=sigma\sqrt{-2log(1-u1)}, theta=2piu2. また，sgima^2=1/2 の場合は r=\sqrt{-log(1-u1)} */
    //sigma = 1/sqrt(2.0);
    //r[i] = sigma * sqrt(-2.0 * log10(1.0-u1[i]));
    r = sigma * sqrt(-2.0 * log(1.0 - u1));
    theta = 2.0 * M_PI * u2;

    //(手順3) x=r*cos*theta, y=r*sin*thetaとする
    *In = r * cos(theta);
    *Qu = r * sin(theta);

    return;

}
