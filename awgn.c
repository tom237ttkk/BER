//�K�E�X�ϐ��̐���
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include  "MT.h"

#define GENRAND_MAX 0xffffffffL //MT�@�ɂ�萶������闐���̍ő�l4294967295
#define gnrnd() ((1.0/(GENRAND_MAX+1.0))*genrand_int32()) // ���������}�N��

void awgn(double, double*, double*);

void awgn(double sigma, double* In, double* Qu) {

    //�ϐ��錾
    double u1, u2, r, theta; //���������p�ϐ�

    //(�菇1) ���[0,1)��̈�l�����̐���
    u1 = gnrnd();
    u2 = gnrnd();

    /*(�菇2) r=sigma\sqrt{-2log(1-u1)}, theta=2piu2. �܂��Csgima^2=1/2 �̏ꍇ�� r=\sqrt{-log(1-u1)} */
    //sigma = 1/sqrt(2.0);
    //r[i] = sigma * sqrt(-2.0 * log10(1.0-u1[i]));
    r = sigma * sqrt(-2.0 * log(1.0 - u1));
    theta = 2.0 * M_PI * u2;

    //(�菇3) x=r*cos*theta, y=r*sin*theta�Ƃ���
    *In = r * cos(theta);
    *Qu = r * sin(theta);

    return;

}
