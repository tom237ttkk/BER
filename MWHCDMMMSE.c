//MWHCDMMMSE�̌v�Z�@�V�~�����[�V����
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include  "MT.h"
#include <limits.h>
#include <float.h>

#define N 1e6 //�T���v����
#define DB 14 //EbN0�͈̔�
#define GENRAND_MAX 0xffffffffL //MT�@�ɂ�萶������闐���̍ő�l4294967295
#define gnrnd() ((1.0/(GENRAND_MAX+1.0))*genrand_int32()) // ���������}�N��

#pragma warning( disable : 4996)

void awgn(double, double*, double*);
void fwht(int, int, int, double, double*, double*, double*, double*);
int BinXor(int, int);

int main(void) {

    //�ϐ��錾
    int i, j, k, loop; //�J�E���^�ϐ�
    unsigned long Nl;
    int n, m, nh;   double rrn; //fwht�p�ϐ�
    int Nt, Nt2; //���d��
    int Nc, Nd, Np; //�����Ւf�p�p�����[�^
    double Eb, Ao, Gam; //DB�p�ϐ�
    double rho, No, sigma, nI, nQ; //AWGN�p�ϐ�
    int NEb, Neb, Nes; //BER�J�E���^�ϐ�
    double Pb, PFth, BPSKth, BPSKth2; //BER�v�Z�p�ϐ�
    int Nmin; //�ŏ��r�b�g��萔

    No = 1.0; //�G���̓d�̓X�y�N�g�����x
    sigma = sqrt(No / 2.0); //�G���̕��U
    
    clock_t s0clock, e0clock;

    unsigned seed;
    seed = time(NULL);
    init_genrand(seed);

    FILE* fps;
    if ((fps = fopen("seed.dat", "w")) == NULL) {
        printf("\a �t�@�C�� seed.dat ���I�[�v���ł��܂���\n");
        exit(1);
    }
    fprintf(fps, "%d\n", seed);
    fclose(fps);

    printf("seed : %d\n", seed);

    s0clock = clock();
    for (j = 0; j < 1; j++) {
        //Tp=38.9msec�ŌŒ�
        //rho=0.086�ł�Td=rhoTp=3.3454msec
        //rho=0.321�ł�Td=rhoTp=12.4869msec
        //j=0,1��Rc=64kbps��j=2,3��Rc=384kpbs�̏ꍇ

        if (j == 0) {
            rho = 0.086;
            Np = 2490; //Rc=64kbps�ł�Np=RcTp=2489.6 ...
            Nd = 215; //Rc=64kbps�ł�Nd=RcTd=214.10 ...
            Nt = 2048; m = 11;
        }
        else if (j == 1) {
            rho = 0.321;
            Np = 2490; //Rc=64kbps�ł�Np=RcTp=2489.5 ...
            Nd = 800; //Rc=64kbps�ł�Nd=RcTd=799.16 ...
            Nt = 2048; m = 11;
        }
        else if (j == 2) {
            rho = 0.086;
            Np = 14938; //Rc=384kbps�ł�Np=RcTp=14937.6 ...
            Nd = 1285; //Rc=384kbps�ł�Nd=RcTd=1284.6 ...
            Nt = 16384; m = 14;
        } 
        else if (j == 3) {
            rho = 0.321;
            Np = 14938; //Rc=384kbps�ł�Np=RcTp=14937.6 ...
            Nd = 4795; //Rc=384kbps�ł�Nd=RcTd=4794.9 ...
            Nt = 16384; m = 14;
        }
        

        printf("���d��N=%d, �Վ���Np=%d, �Ւf��Nd=%d, �Ւf��=%lf, N/Np=%lf\n", Nt, Np, Nd, rho, (double)Nt / (double)Np);
        printf("���r�b�g��N_EB / ���V���{����N_ES / dB�l / BER / ���EBER / ���d���Ȃ�BER / BPSK���_�l\n");

        Nt2 = Nt / 2;
        n = (int)pow(2.0, (double)m);
        nh = n / 2;
        rrn = 1.0 / sqrt((double)n);

        int* a;//���M���p�ϐ�
        a = calloc(Nt, sizeof(double));
        if (a == NULL) {
            puts("�z��a�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        double* bRe, * bIm; //�ꎟ�ϒ��M��
        bRe = calloc(Nt, sizeof(double));  bIm = calloc(Nt, sizeof(double));
        if (bRe == NULL || bIm == NULL) {
            puts("�z��bRe��bIm�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        double* sRe, * sIm; //WHCDM���d���㑗�M�M��
        sRe = calloc(Nt, sizeof(double));  sIm = calloc(Nt, sizeof(double));
        if (sRe == NULL || sIm == NULL) {
            puts("�z��sRe��sIm�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        double* hRe, * hIm; //�`���l�����p�ϐ�
        hRe = calloc(Nt, sizeof(double));  hIm = calloc(Nt, sizeof(double));
        if (hRe == NULL || hIm == NULL) {
            puts("�z��hRe��hIm�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        double* rRe, * rIm; //��M�M���p�ϐ�
        rRe = calloc(Nt, sizeof(double));  rIm = calloc(Nt, sizeof(double));
        if (rRe == NULL || rIm == NULL) {
            puts("�z��rRe��rIm�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        double* rcRe, * rcIm; //������M�M���p�ϐ�
        rcRe = calloc(Nt2, sizeof(double));  rcIm = calloc(Nt2, sizeof(double));
        if (rcRe == NULL || rcIm == NULL) {
            puts("�z��rcRe��rcIm�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        double* wRe, * wIm; //�d�݂Â��W��
        wRe = calloc(Nt, sizeof(double));   wIm = calloc(Nt, sizeof(double));
        if (wRe == NULL || wIm == NULL) {
            puts("�z��wRe��wIm�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        double* dRe, * dIm; //���d������p�ϐ�
        dRe = calloc(Nt2, sizeof(double));  dIm = calloc(Nt2, sizeof(double));
        if (dRe == NULL || dIm == NULL) {
            puts("�z��dRe��dIm�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }

        int* ah;//�����M���p�ϐ�
        ah = calloc(Nt, sizeof(double));
        if (ah == NULL) {
            puts("�z��ah�̋L���̈�̊m�ۂɎ��s");
            exit(1);
        }


        for (k = 0; k < DB; k++) {
            Eb = pow(10.0, (double)k / 10.0); //�M���̓d��
            Ao = sqrt(Eb / No); //�M���̐U��
            Gam = No / Eb;
            Nes = 0;
            NEb = 0;

            if ( k < 10 ) {
                Nmin = 2000;
            }
            else if ((10 <= k) && (k < 11)) {
                Nmin = 500;
            }
            else if ((11 <= k) && (k < 12)) {
                Nmin = 200;
            }
            else if ((12 <= k) && (k < 13)) {
                Nmin = 100;
            }
            else if ( 13 <= k ) {
                Nmin = 10;
            }

            for (Nl = 0; Nes < Nmin; Nl++) {
                Nc = 0;     Neb = 0;
                for (loop = 0; loop < Np; loop++) {
                    //��񌹐���
                    for (i = 0; i < Nt2; i++) {
                        //1.�f�[�^�����@�\
                        //k = 2.0 * gnrnd() �Ƃ���
                        //k<1�Ȃ�a[i]=0�ŏ��0�����1�Ƃ��ă}�b�s���O
                        //1<k�Ȃ�a[i]=1�ŏ��1�����-1�Ƃ��ă}�b�s���O
                        a[i] = (int)(2.0 * gnrnd());
                        a[Nt2 + i] = (int)(2.0 * gnrnd());

                        //2.�ꎟ�ϒ��@�\
                        //�O�����̓��A���̂�
                        bRe[i] = Ao * cos(((double)a[i]) * M_PI);
                        bIm[i] = 0.0;
                        //�㔼���̓C���[�W�̂�
                        bRe[Nt2 + i] = 0.0;
                        bIm[Nt2 + i] = Ao * cos(((double)a[Nt2 + i]) * M_PI);
                    }
                    //����ɂ�� b = bRe + j*bIm �ŕ��f�M���𐶐������̂Ɠ���
                    //
                    /* 0<=i<Nt2�ɂ�����
                    bRe[i]=0�ȊO, bIm[Nt2+i]=0, bIm[i]=0, bIm[Nt2+i]=0�ȊO������
                    b[i] = bRe[i] + j*bIm[i] = bRe[i],
                    b[Nt2+i] = bRe[Nt2+i] + j*bIm[Nt2+i] = j*bIm[Nt2+i] */

                    //3.���d���@�\
                    fwht(m, n, nh, rrn, bRe, bIm, sRe, sIm);
                    //free(bRe); free(bIm);//�������J��
                    //MWHCDM�ɂ�� bRe to sRe,   bIm to sIm (WHCDM�̌W���͎����݂̂̂��ߎ����Ƌ����𕪂��Čv�Z�\)
                    //�O������ s[i] = sRe[i] - j*sIm[i] = bRe[i] + j+bIm[Nt2+i] (�W���͏ȗ�)
                    //�㔼���� s[Nt2+i] = sRe[Nt2+i] - j*sIm[Nt2+i] = bRe[i] - j*bIm[Nt2+i] (�W���͏ȗ�)

                    //4.�����I�Ւf��AWGN���󂯂���M�M��
                    for (i = 0; i < Nt; i++) {
                        // Nc�̓`�b�v�̈ʒu��\��
                        if (Nc < Nd) { //�M�����Ւf�����ꍇ
                            hRe[i] = 0.0; //�Ւf�̓`���l������0
                            hIm[i] = 0.0;
                        }
                        else { //�M�����Ւf����Ȃ��ꍇ
                            hRe[i] = 1.0; //�Ւf�Ȃ����̓`���l������1
                            hIm[i] = 0.0; //�`���l���̕��f�v�f�͍���0�ōl����
                        }
                        awgn(sigma, &nI, &nQ); //AWGN
                        rRe[i] = (hRe[i] * sRe[i] - hIm[i] * sIm[i]) + nI;
                        rIm[i] = (hRe[i] * sIm[i] + hIm[i] * sRe[i]) + nQ;
                        Nc = (Nc + 1) % Np;
                    }
                    //free(sRe); free(sIm);//�������J��

                    //5.������M�M���̍쐬��MMSE�d�݂Â�              
                    for (i = 0; i < Nt2; i++) {
                        double Den;
                        Den = (pow(hRe[i], 2.0) + pow(hIm[i], 2.0)) + (pow(hRe[Nt2 + i], 2.0) + pow(hIm[Nt2 + i], 2.0)) + Gam;
                        //�O�����W��
                        wRe[i] = hRe[i] / Den;
                        wIm[i] = -1.0 * hIm[i] / Den;
                        //�㔼���W��
                        wRe[Nt2 + i] = hRe[Nt2 + i] / Den;
                        wIm[Nt2 + i] = hIm[Nt2 + i] / Den;
                        //������M�M������
                        rIm[Nt2 + i] = -1.0 * rIm[Nt2 + i]; //�㔼�̎�M�M���̂݋������Ƃ�
                        rcRe[i] = (wRe[i] * rRe[i] - wIm[i] * rIm[i]) + (wRe[Nt2 + i] * rRe[Nt2 + i] - wIm[Nt2 + i] * rIm[Nt2 + i]);
                        rcIm[i] = (wRe[i] * rIm[i] + wIm[i] * rRe[i]) + (wRe[Nt2 + i] * rIm[Nt2 + i] + wIm[Nt2 + i] * rRe[Nt2 + i]);
                    }
                    //free(rRe); free(rIm); free(hRe); free(hIm); free(wRe); free(wIm); //�������J��
                    // rc = w[i]*r[i] + w[i]*r[i]^{*} 
                    //�M������ Nt2=Nt/2

                    //6.���d�����@�\
                    double n2, nh2, rrn2;
                    n2 = (int)pow(2.0, (double)(m - 1));
                    nh2 = n2 / 2.0;
                    rrn2 = 1.0 / sqrt((double)n);
                    fwht((m - 1), n2, nh2, rrn2, rcRe, rcIm, dRe, dIm);
                    //free(rcRe); free(rcIm); //�������J��
                    //�ŏI�I�� D = [ d0  d1 ] = [ dRe  dIm ]

                    //7.�f�[�^����E�r�b�g��萔�J�E���g�@�\
                    for (i = 0; i < Nt2; i++) {
                        if (dRe[i] >= 0) { ah[i] = 0; }
                        else { ah[i] = 1; }
                        if (dIm[i] >= 0) { ah[Nt2 + i] = 0; }
                        else { ah[Nt2 + i] = 1; }
                        Neb = Neb + BinXor(a[i], ah[i]) + BinXor(a[Nt2 + i], ah[Nt2 + i]);
                    }
                    //free(dRe); free(dIm); //�������J��
                    NEb = NEb + Neb;
                    if (Neb > 0) Nes = Nes + 1;
                    Neb = 0;
                }
            }


            //BER�v�Z
            Pb = (double)NEb / ((double)Np * (double)Nl * (double)Nt);
            PFth = (1.0 / 2.0) * (rho + (1.0 - rho) * erfc(sqrt(Eb / No)));
            BPSKth = (1.0 / 2.0) * erfc(sqrt(Eb / No));
            BPSKth2 = (1.0 / 2.0) * erfc(sqrt((1.0 - rho) * Eb / No));

            FILE* fp0, * fp1, * fp2, * fp3;
            if ((fp0 = fopen("MWHCDMMMSE384kbps_0086.dat", "a")) == NULL) {
                printf("\a �t�@�C�� MWHCDM.dat ���I�[�v���ł��܂���\n");
                exit(1);
            }
            if ((fp1 = fopen("woMutiplexing_rho0086.dat", "a")) == NULL) {
                printf("\a �t�@�C�� woMultiplexing.dat ���I�[�v���ł��܂���\n");
                exit(1);
            }
            if ((fp2 = fopen("BPSKth.dat", "a")) == NULL) {
                printf("\a �t�@�C�� BPSKth.dat ���I�[�v���ł��܂���\n");
                exit(1);
            }
            if ((fp3 = fopen("BPSKth_rho0086.dat", "a")) == NULL) {
                printf("\a �t�@�C�� BPSKth_rho.dat ���I�[�v���ł��܂���\n");
                exit(1);
            }

            fprintf(fp0, "%d \t %8.3e\n", k, Pb);
            fprintf(fp1, "%d \t %8.3e\n", k, PFth);
            fprintf(fp2, "%d \t %8.3e\n", k, BPSKth);
            fprintf(fp3, "%d \t %8.3e\n", k, BPSKth2);

            fclose(fp0); fclose(fp1); fclose(fp2); fclose(fp3);

            printf("%d / %d / %d / %8.3e / %8.3e / %8.3e / %8.3e\n", NEb, Nes, k, Pb, BPSKth2, PFth, BPSKth);

        }
        free(a); free(bRe); free(bIm); free(sRe); free(sIm); free(rRe); free(rIm); free(hRe); free(hIm); free(wRe); free(wIm);
        free(rcRe); free(rcIm); free(dRe); free(dIm); free(ah);
    }
    e0clock = clock();
    printf("���s���x:%f\n", ((double)e0clock - (double)s0clock) / CLOCKS_PER_SEC);

    return 0;
}

int BinXor(int v, int w) {
    int h;
    h = abs(v - w);
    return h;
}
