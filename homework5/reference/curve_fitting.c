#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define  M  6    /* データのペア数 */
#define  N  3    /* N次式で近似 */

/* ベクトルの入力 */
void input_vector2( double *b, char c, FILE *fin, FILE *fout);
/*  行列の領域確保 */
/* 部分ピボット選択付きガウス消去法 */
void gauss2( double a[N+1][N+1], double b[N+1], int n ); 
/* 最小2乗近似呼び出し */
void least_square_coef( double x[M], double y[M], double a[N+1],  FILE *fout );
/* 最小2乗近似 */
void least_square( double x[M], double y[M], FILE *fout );
/* 曲線の生成 */
void curve( double x[M], double y[M], double t[M], FILE *fout);
/* 曲線上の点列データ出力 */
void curve_shape( double a0[N+1], double b0[N+1], FILE *fp1);

double a[N+1], b[N+1];

int main(void)
{
  FILE *fin, *fout, *fout1;
  double x[M], y[M], t[M]; 
  errno_t error;

  /* ファイルのオープン */
  if ( (error = fopen_s(&fin, "curve_fitting.dat", "r")) != 0 )
  {
      printf("ファイルが見つかりません : curve_fitting.dat \n");
      exit(1);
  }
  if( (error = fopen_s(&fout, "output_curve.dat", "w")) != 0 )
  {
    printf("ファイルが作成できません : output_curve.dat \n");
    exit(1);
   }
  if( (error = fopen_s(&fout1, "output_shape.dat", "w")) != 0 )
  {
    printf("ファイルが作成できません : output_shape.dat \n");
    exit(1);
   }

  input_vector2( x, 'x', fin, fout );    /* ベクトルxの入出力 */
  input_vector2( y, 'y', fin, fout );    /* ベクトルyの入出力 */
  input_vector2( t, 't', fin, fout );    /* パラメータtの入出力 */

  curve( x, y, t, fout);  /* 最小2乗近似 */

  curve_shape( a, b, fout1);

  fclose(fin); fclose(fout);  /* ファイルのクローズ */

  return 0;
}

void curve( double x[M], double y[M], double t[M], FILE *fout)
{
	least_square_coef( t, x, a, fout);
	least_square_coef( t, y, b, fout);
}

void least_square_coef( double x[M], double y[M], double a[N+1],  FILE *fout )
{
  double p[N+1][N+1];
  int i, j, k;

  /* 右辺ベクトルの作成 */
  for(i=0; i <= N; i++)
  {
    a[i]=0.0;
    for( j = 0; j < M; j++)
    {
      a[i] +=  y[j]*pow(x[j],(double)(i)) ;
    }
  }

  /* 係数行列の作成 */
  for( i = 0; i <= N; i++)
  {
    for( j = 0; j <= i;  j++ )
    {
      p[i][j]=0.0;
      for( k =0; k < M; k++)
      {
        p[i][j] += pow( x[k], (double)(i+j) );
      }
      p[j][i] = p[i][j];
    } 
}
  /* 連立1次方程式を解く. 結果はaに上書き */
   gauss2( p, a, N+1 );            

  /* 結果の出力 */
  fprintf( fout, "最小2乗近似式はy=");
  for( i = N ; i >= 0 ; i--)
  {
		if(i==N){
			fprintf(fout, "%5.2f x^%d ", a[i],i);
		}
		else{
			if(a[i]>0){
				fprintf(fout, "+ %5.2f x^%d ", a[i],i);
			}
			else{
				fprintf(fout, "- %5.2f x^%d ",fabs(a[i]),i);
			}
		}
  }
  fprintf(fout, "\n"); 
}

/* 部分ピボット選択付きガウス消去法 */
void gauss2( double a[N+1][N+1], double b[N+1], int n )
{
  int i, j, k, ip;
  double alpha, tmp;
  double amax, eps=pow(2.0, -50.0); /* eps = 2^{-50}とする */

  for( k = 0; k < n-1; k++)
  {
    /* ピボットの選択 */
    amax = fabs(a[k][k]); ip = k;
    for( i = k+1; i < n; i++)
    {
      if ( fabs(a[i][k]) > amax )
      {
        amax = fabs(a[i][k]); ip = i;
      }
    }
    /* 正則性の判定 */
    if ( amax < eps ) printf("入力した行列は正則ではない!!\n");
    /* 行交換 */
    if ( ip != k)
    {
      for( j = k; j < n; j++)
      {
        tmp = a[k][j]; a[k][j]=a[ip][j]; a[ip][j]=tmp;
      }
        tmp = b[k] ; b[k]=b[ip]; b[ip]=tmp;
    }
      /* 前進消去 */
    for( i = k+1; i < n; i++)
    {
      alpha = - a[i][k]/a[k][k];
      for( j = k+1; j < n; j++)
      {
        a[i][j] = a[i][j] + alpha * a[k][j];
      }
      b[i] = b[i] + alpha * b[k];
    }
  }

  /* 後退代入 */
  b[n-1] = b[n-1]/a[n-1][n-1];
  for( k = n-2; k >= 0; k--)
  {
    tmp = b[k];
    for( j = k+1; j <= n; j++)
    {
      tmp = tmp - a[k][j] * b[j];
    }
    b[k] = tmp/a[k][k];
  }
}

/* b[1...n]の入力 */
void input_vector2( double b[M], char c, FILE *fin, FILE *fout)
{
  int i;

  fprintf( fout, "ベクトル%cは次の通りです\n", c);
  for( i = 0 ; i < M ; i++)
  {
    fscanf_s(fin, "%lf", &b[i]);
    fprintf(fout, "%5.2f\t", b[i]);
  }
  fprintf( fout, "\n");
}

void curve_shape( double a0[N+1], double b0[N+1], FILE *fp1)
{
	int i,j;
	double dt0=0.05;
	double x0,y0,t0=0.,t1;
	for(i=0;i<=20;++i){
		x0=a0[0];
		y0=b0[0];
		t0+=dt0;
		t1=t0;
		for(j=1;j<=N;++j){
			x0+=a0[j]*t1;
			y0+=b0[j]*t1;
			t1*=t0;
		}
		fprintf(fp1,"%lf\t%lf\n",x0,y0);
	}
}
