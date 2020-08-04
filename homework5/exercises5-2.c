#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//cl /source-charset:utf-8 exercises5-2.c /EHsc

#define M 11 /* データのペア数 */
#define N 5  /* N次式で近似 */

/* ベクトルの入力 */
void input_vector2(double *b, char c, int n, FILE *fin, FILE *fout);
/* 部分ピボット選択付きガウス消去法 */
void gauss2(double a[N + 1][N + 1], double b[N + 1], int n);
/* 最小2乗近似 */
void least_square(double *x, double *y, FILE *fout);

int main(void)
{
    FILE *fin, *fout;
    double x[M], y[M];
    errno_t error;

    /* ファイルのオープン */
    if ((error = fopen_s(&fin, "input_func5-1.dat", "r")) != 0)
    {
        printf("ファイルが見つかりません : input_func.dat \n");
        exit(1);
    }
    if ((error = fopen_s(&fout, "output_func5-2-5zi.dat", "w")) != 0)
    {
        printf("ファイルが作成できません : output_func.dat \n");
        exit(1);
    }

    input_vector2(x, 'x', M, fin, fout); /* ベクトルxの入出力 */
    input_vector2(y, 'y', M, fin, fout); /* ベクトルyの入出力 */

    least_square(x, y, fout); /* 最小2乗近似 */

    fclose(fin);
    fclose(fout); /* ファイルのクローズ */

    return 0;
}

void least_square(double x[M], double y[M], FILE *fout)
{
    double a[N + 1], p[N + 1][N + 1];
    int i, j, k;

    /* 右辺ベクトルの作成 */
    for (i = 0; i <= N; i++)
    {
        a[i] = 0.0;
        for (j = 0; j < M; j++)
        {
            a[i] += y[j] * pow(x[j], (double)i);
        }
    }

    /* 係数行列の作成 */
    for (i = 0; i <= N; i++)
    {
        for (j = 0; j <= i; j++)
        {
            p[i][j] = 0.0;
            for (k = 0; k < M; k++)
            {
                p[i][j] += pow(x[k], (double)(i + j));
            }
            p[j][i] = p[i][j];
        }
    }
    /* 連立1次方程式を解く. 結果はaに上書き */
    gauss2(p, a, N + 1);

    /* 結果の出力 */
    fprintf(fout, "最小2乗近似式は y =");
    for (i = N; i >= 0; i--)
    {
        if (a[i] > 0)
        {
            if (i == N)
            {
                fprintf(fout, " %5.2f x^%d ", a[i], i);
            }
            else
            {
                fprintf(fout, "+ %5.2f x^%d ", a[i], i);
            }
        }
        else
        {
            fprintf(fout, "- %5.2f x^%d ", fabs(a[i]), i);
        }
    }
    fprintf(fout, "\n");
}

/* 部分ピボット選択付きガウス消去法 */
void gauss2(double a[N + 1][N + 1], double b[N + 1], int n)
{
    int i, j, k, ip;
    double alpha, tmp;
    double amax, eps = pow(2.0, -50.0); /* eps = 2^{-50}とする */

    for (k = 0; k < n - 1; k++)
    {
        /* ピボットの選択 */
        amax = fabs(a[k][k]);
        ip = k;
        for (i = k + 1; i < n; i++)
        {
            if (fabs(a[i][k]) > amax)
            {
                amax = fabs(a[i][k]);
                ip = i;
            }
        }
        /* 正則性の判定 */
        if (amax < eps)
            printf("入力した行列は正則ではない!!\n");
        /* 行交換 */
        if (ip != k)
        {
            for (j = k; j < n; j++)
            {
                tmp = a[k][j];
                a[k][j] = a[ip][j];
                a[ip][j] = tmp;
            }
            tmp = b[k];
            b[k] = b[ip];
            b[ip] = tmp;
        }
        /* 前進消去 */
        for (i = k + 1; i < n; i++)
        {
            alpha = -a[i][k] / a[k][k];
            for (j = k + 1; j < n; j++)
            {
                a[i][j] = a[i][j] + alpha * a[k][j];
            }
            b[i] = b[i] + alpha * b[k];
        }
    }

    /* 後退代入 */
    b[n - 1] = b[n - 1] / a[n - 1][n - 1];
    for (k = n - 2; k >= 0; k--)
    {
        tmp = b[k];
        for (j = k + 1; j < n; j++)
        {
            tmp = tmp - a[k][j] * b[j];
        }
        b[k] = tmp / a[k][k];
    }
}

/* b[1...n]の入力 */
void input_vector2(double b[N + 1], char c, int n, FILE *fin, FILE *fout)
{
    int i;

    fprintf(fout, "ベクトル%cは次の通りです\n", c);
    for (i = 0; i < n; i++)
    {
        fscanf_s(fin, "%lf", &b[i]);
        fprintf(fout, "%5.2f\t", b[i]);
    }
    fprintf(fout, "\n");
}
