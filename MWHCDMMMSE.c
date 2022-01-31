//MWHCDMMMSEの計算機シミュレーション
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include  "MT.h"
#include <limits.h>
#include <float.h>

#define N 1e6 //サンプル数
#define DB 14 //EbN0の範囲
#define GENRAND_MAX 0xffffffffL //MT法により生成される乱数の最大値4294967295
#define gnrnd() ((1.0/(GENRAND_MAX+1.0))*genrand_int32()) // 乱数生成マクロ

#pragma warning( disable : 4996)

void awgn(double, double*, double*);
void fwht(int, int, int, double, double*, double*, double*, double*);
int BinXor(int, int);

int main(void) {

    //変数宣言
    int i, j, k, loop; //カウンタ変数
    unsigned long Nl;
    int n, m, nh;   double rrn; //fwht用変数
    int Nt, Nt2; //多重数
    int Nc, Nd, Np; //周期遮断用パラメータ
    double Eb, Ao, Gam; //DB用変数
    double rho, No, sigma, nI, nQ; //AWGN用変数
    int NEb, Neb, Nes; //BERカウンタ変数
    double Pb, PFth, BPSKth, BPSKth2; //BER計算用変数
    int Nmin; //最小ビット誤り数

    No = 1.0; //雑音の電力スペクトル密度
    sigma = sqrt(No / 2.0); //雑音の分散
    
    clock_t s0clock, e0clock;

    unsigned seed;
    seed = time(NULL);
    init_genrand(seed);

    FILE* fps;
    if ((fps = fopen("seed.dat", "w")) == NULL) {
        printf("\a ファイル seed.dat をオープンできません\n");
        exit(1);
    }
    fprintf(fps, "%d\n", seed);
    fclose(fps);

    printf("seed : %d\n", seed);

    s0clock = clock();
    for (j = 0; j < 1; j++) {
        //Tp=38.9msecで固定
        //rho=0.086ではTd=rhoTp=3.3454msec
        //rho=0.321ではTd=rhoTp=12.4869msec
        //j=0,1はRc=64kbpsのj=2,3はRc=384kpbsの場合

        if (j == 0) {
            rho = 0.086;
            Np = 2490; //Rc=64kbpsではNp=RcTp=2489.6 ...
            Nd = 215; //Rc=64kbpsではNd=RcTd=214.10 ...
            Nt = 2048; m = 11;
        }
        else if (j == 1) {
            rho = 0.321;
            Np = 2490; //Rc=64kbpsではNp=RcTp=2489.5 ...
            Nd = 800; //Rc=64kbpsではNd=RcTd=799.16 ...
            Nt = 2048; m = 11;
        }
        else if (j == 2) {
            rho = 0.086;
            Np = 14938; //Rc=384kbpsではNp=RcTp=14937.6 ...
            Nd = 1285; //Rc=384kbpsではNd=RcTd=1284.6 ...
            Nt = 16384; m = 14;
        } 
        else if (j == 3) {
            rho = 0.321;
            Np = 14938; //Rc=384kbpsではNp=RcTp=14937.6 ...
            Nd = 4795; //Rc=384kbpsではNd=RcTd=4794.9 ...
            Nt = 16384; m = 14;
        }
        

        printf("多重数N=%d, 遮周期Np=%d, 遮断長Nd=%d, 遮断率=%lf, N/Np=%lf\n", Nt, Np, Nd, rho, (double)Nt / (double)Np);
        printf("誤りビット数N_EB / 誤りシンボル数N_ES / dB値 / BER / 下界BER / 多重化なしBER / BPSK理論値\n");

        Nt2 = Nt / 2;
        n = (int)pow(2.0, (double)m);
        nh = n / 2;
        rrn = 1.0 / sqrt((double)n);

        int* a;//源信号用変数
        a = calloc(Nt, sizeof(double));
        if (a == NULL) {
            puts("配列aの記憶領域の確保に失敗");
            exit(1);
        }

        double* bRe, * bIm; //一次変調信号
        bRe = calloc(Nt, sizeof(double));  bIm = calloc(Nt, sizeof(double));
        if (bRe == NULL || bIm == NULL) {
            puts("配列bReとbImの記憶領域の確保に失敗");
            exit(1);
        }

        double* sRe, * sIm; //WHCDM多重化後送信信号
        sRe = calloc(Nt, sizeof(double));  sIm = calloc(Nt, sizeof(double));
        if (sRe == NULL || sIm == NULL) {
            puts("配列sReとsImの記憶領域の確保に失敗");
            exit(1);
        }

        double* hRe, * hIm; //チャネル情報用変数
        hRe = calloc(Nt, sizeof(double));  hIm = calloc(Nt, sizeof(double));
        if (hRe == NULL || hIm == NULL) {
            puts("配列hReとhImの記憶領域の確保に失敗");
            exit(1);
        }

        double* rRe, * rIm; //受信信号用変数
        rRe = calloc(Nt, sizeof(double));  rIm = calloc(Nt, sizeof(double));
        if (rRe == NULL || rIm == NULL) {
            puts("配列rReとrImの記憶領域の確保に失敗");
            exit(1);
        }

        double* rcRe, * rcIm; //合成受信信号用変数
        rcRe = calloc(Nt2, sizeof(double));  rcIm = calloc(Nt2, sizeof(double));
        if (rcRe == NULL || rcIm == NULL) {
            puts("配列rcReとrcImの記憶領域の確保に失敗");
            exit(1);
        }

        double* wRe, * wIm; //重みづけ係数
        wRe = calloc(Nt, sizeof(double));   wIm = calloc(Nt, sizeof(double));
        if (wRe == NULL || wIm == NULL) {
            puts("配列wReとwImの記憶領域の確保に失敗");
            exit(1);
        }

        double* dRe, * dIm; //多重分離後用変数
        dRe = calloc(Nt2, sizeof(double));  dIm = calloc(Nt2, sizeof(double));
        if (dRe == NULL || dIm == NULL) {
            puts("配列dReとdImの記憶領域の確保に失敗");
            exit(1);
        }

        int* ah;//復調信号用変数
        ah = calloc(Nt, sizeof(double));
        if (ah == NULL) {
            puts("配列ahの記憶領域の確保に失敗");
            exit(1);
        }


        for (k = 0; k < DB; k++) {
            Eb = pow(10.0, (double)k / 10.0); //信号の電力
            Ao = sqrt(Eb / No); //信号の振幅
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
                    //情報源生成
                    for (i = 0; i < Nt2; i++) {
                        //1.データ生成機能
                        //k = 2.0 * gnrnd() として
                        //k<1ならa[i]=0で情報源0を後に1としてマッピング
                        //1<kならa[i]=1で情報源1を後に-1としてマッピング
                        a[i] = (int)(2.0 * gnrnd());
                        a[Nt2 + i] = (int)(2.0 * gnrnd());

                        //2.一次変調機能
                        //前半部はリアルのみ
                        bRe[i] = Ao * cos(((double)a[i]) * M_PI);
                        bIm[i] = 0.0;
                        //後半部はイメージのみ
                        bRe[Nt2 + i] = 0.0;
                        bIm[Nt2 + i] = Ao * cos(((double)a[Nt2 + i]) * M_PI);
                    }
                    //これにより b = bRe + j*bIm で複素信号を生成したのと等価
                    //
                    /* 0<=i<Nt2において
                    bRe[i]=0以外, bIm[Nt2+i]=0, bIm[i]=0, bIm[Nt2+i]=0以外だから
                    b[i] = bRe[i] + j*bIm[i] = bRe[i],
                    b[Nt2+i] = bRe[Nt2+i] + j*bIm[Nt2+i] = j*bIm[Nt2+i] */

                    //3.多重化機能
                    fwht(m, n, nh, rrn, bRe, bIm, sRe, sIm);
                    //free(bRe); free(bIm);//メモリ開放
                    //MWHCDMにより bRe to sRe,   bIm to sIm (WHCDMの係数は実数のみのため実部と虚部を分けて計算可能)
                    //前半部は s[i] = sRe[i] - j*sIm[i] = bRe[i] + j+bIm[Nt2+i] (係数は省略)
                    //後半部は s[Nt2+i] = sRe[Nt2+i] - j*sIm[Nt2+i] = bRe[i] - j*bIm[Nt2+i] (係数は省略)

                    //4.周期的遮断とAWGNを受けた受信信号
                    for (i = 0; i < Nt; i++) {
                        // Ncはチップの位置を表す
                        if (Nc < Nd) { //信号が遮断される場合
                            hRe[i] = 0.0; //遮断はチャネル利得0
                            hIm[i] = 0.0;
                        }
                        else { //信号が遮断されない場合
                            hRe[i] = 1.0; //遮断なし時はチャネル利得1
                            hIm[i] = 0.0; //チャネルの複素要素は今は0で考える
                        }
                        awgn(sigma, &nI, &nQ); //AWGN
                        rRe[i] = (hRe[i] * sRe[i] - hIm[i] * sIm[i]) + nI;
                        rIm[i] = (hRe[i] * sIm[i] + hIm[i] * sRe[i]) + nQ;
                        Nc = (Nc + 1) % Np;
                    }
                    //free(sRe); free(sIm);//メモリ開放

                    //5.合成受信信号の作成とMMSE重みづけ              
                    for (i = 0; i < Nt2; i++) {
                        double Den;
                        Den = (pow(hRe[i], 2.0) + pow(hIm[i], 2.0)) + (pow(hRe[Nt2 + i], 2.0) + pow(hIm[Nt2 + i], 2.0)) + Gam;
                        //前半部係数
                        wRe[i] = hRe[i] / Den;
                        wIm[i] = -1.0 * hIm[i] / Den;
                        //後半部係数
                        wRe[Nt2 + i] = hRe[Nt2 + i] / Den;
                        wIm[Nt2 + i] = hIm[Nt2 + i] / Den;
                        //合成受信信号生成
                        rIm[Nt2 + i] = -1.0 * rIm[Nt2 + i]; //後半の受信信号のみ共役をとる
                        rcRe[i] = (wRe[i] * rRe[i] - wIm[i] * rIm[i]) + (wRe[Nt2 + i] * rRe[Nt2 + i] - wIm[Nt2 + i] * rIm[Nt2 + i]);
                        rcIm[i] = (wRe[i] * rIm[i] + wIm[i] * rRe[i]) + (wRe[Nt2 + i] * rIm[Nt2 + i] + wIm[Nt2 + i] * rRe[Nt2 + i]);
                    }
                    //free(rRe); free(rIm); free(hRe); free(hIm); free(wRe); free(wIm); //メモリ開放
                    // rc = w[i]*r[i] + w[i]*r[i]^{*} 
                    //信号長は Nt2=Nt/2

                    //6.多重分離機能
                    double n2, nh2, rrn2;
                    n2 = (int)pow(2.0, (double)(m - 1));
                    nh2 = n2 / 2.0;
                    rrn2 = 1.0 / sqrt((double)n);
                    fwht((m - 1), n2, nh2, rrn2, rcRe, rcIm, dRe, dIm);
                    //free(rcRe); free(rcIm); //メモリ開放
                    //最終的に D = [ d0  d1 ] = [ dRe  dIm ]

                    //7.データ判定・ビット誤り数カウント機能
                    for (i = 0; i < Nt2; i++) {
                        if (dRe[i] >= 0) { ah[i] = 0; }
                        else { ah[i] = 1; }
                        if (dIm[i] >= 0) { ah[Nt2 + i] = 0; }
                        else { ah[Nt2 + i] = 1; }
                        Neb = Neb + BinXor(a[i], ah[i]) + BinXor(a[Nt2 + i], ah[Nt2 + i]);
                    }
                    //free(dRe); free(dIm); //メモリ開放
                    NEb = NEb + Neb;
                    if (Neb > 0) Nes = Nes + 1;
                    Neb = 0;
                }
            }


            //BER計算
            Pb = (double)NEb / ((double)Np * (double)Nl * (double)Nt);
            PFth = (1.0 / 2.0) * (rho + (1.0 - rho) * erfc(sqrt(Eb / No)));
            BPSKth = (1.0 / 2.0) * erfc(sqrt(Eb / No));
            BPSKth2 = (1.0 / 2.0) * erfc(sqrt((1.0 - rho) * Eb / No));

            FILE* fp0, * fp1, * fp2, * fp3;
            if ((fp0 = fopen("MWHCDMMMSE384kbps_0086.dat", "a")) == NULL) {
                printf("\a ファイル MWHCDM.dat をオープンできません\n");
                exit(1);
            }
            if ((fp1 = fopen("woMutiplexing_rho0086.dat", "a")) == NULL) {
                printf("\a ファイル woMultiplexing.dat をオープンできません\n");
                exit(1);
            }
            if ((fp2 = fopen("BPSKth.dat", "a")) == NULL) {
                printf("\a ファイル BPSKth.dat をオープンできません\n");
                exit(1);
            }
            if ((fp3 = fopen("BPSKth_rho0086.dat", "a")) == NULL) {
                printf("\a ファイル BPSKth_rho.dat をオープンできません\n");
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
    printf("実行速度:%f\n", ((double)e0clock - (double)s0clock) / CLOCKS_PER_SEC);

    return 0;
}

int BinXor(int v, int w) {
    int h;
    h = abs(v - w);
    return h;
}
