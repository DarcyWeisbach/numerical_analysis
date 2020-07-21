#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
//cl /source-charset:shift-jis exercises3-1.c /EHsc
/* 関数の定義 */
double func1(double x);
double func2(double x);
/* 台形公式 */
double trapezoidal( double a, double b, int n, double (*f)(double) );

int main(void)
{
  int n=100;
  double func1_ans=1;
  double func2_ans=M_PI;

  printf("2.0/(x*x)を[1,2]で積分します。分割数は%dです\n", n);
  printf("結果は%20.15fです\n",trapezoidal(1.0, 2.0, n, func1) );
  printf("絶対値誤差は%20.15fです\n",fabs(trapezoidal(1.0,2.0,n,func1)-func1_ans));

  printf("4.0/(1+x*x)を[0,1]で積分します。分割数は%dです\n", n);
  printf("結果は%20.15fです\n",trapezoidal(0.0, 1.0, n, func2) );
  printf("絶対値誤差は%20.15fです\n",fabs(trapezoidal(0.0,1.0,n,func2)-func2_ans));

  return 0;
}

/* 台形公式 */
double trapezoidal( double a, double b, int n, double (*f)(double) )
{
  /*高階関数
    直接関数を呼び出すよりも自由度が高くなる
    */
  double T, h;
  int i;

  h = ( b - a ) /n ;  /* 刻み幅の指定 */

  /* 台形公式 */
  T = ( (*f)(a) + (*f)(b) ) / 2.0;
  for ( i = 1; i < n; i++) T += (*f)( a + i*h );
  T *= h;

  return T;
}

/* 関数の定義 */
double func1(double x)
{
  return( 2.0/(x*x) );
}

double func2(double x)
{
  return( 4.0 / (1.0+x*x) );
}
