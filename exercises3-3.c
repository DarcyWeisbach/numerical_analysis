#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
//cl /source-charset:shift-jis exercises3-3.c /EHsc
/* 関数の定義 */
double func(double x);
/* 台形公式 */
double trapezoidal( double a, double b, int n, double (*f)(double) );
/* シンプソン公式 */
double simpson( double a, double b, int n, double (*f)(double) );

int main(void)
{
  int num_split;
  double func_ans=M_PI;

  printf("4.0*sqrt(1-x*x)を[0,1]で積分します\n");
  printf("\t\t      台形公式(誤差)\t     Simpson公式(誤差)\n");
  for(num_split=2;num_split<=128;num_split*=2){
    printf("分割数%-8d",num_split );
    printf("%5.10f(%-5.6f)",trapezoidal(0.0,1.0,num_split,func),fabs(trapezoidal(0.0,1.0,num_split,func)-func_ans));
    printf("\t%5.10f(%-5.6f)\n",simpson(0.0,1.0,num_split,func),fabs(simpson(0.0,1.0,num_split,func)-func_ans));

  }



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

double simpson( double a, double b, int n, double (*f)(double) )
{
  double S, h;
  int i;
  n=n/(2.0);

  h = ( b - a ) / (2.0*n) ;  /* 刻み幅の指定 */

  /* シンプソン公式 */
  S = ( (*f)(a) + (*f)(b) ) ;
  for ( i = 1; i < n; i++)
  {
    S += 4.0*(*f)( a + (2.0*i-1.0)*h ) + 2.0*(*f)( a + 2.0*i*h );
  }
  S += 4.0*(*f)( a + (2.0*n-1.0)*h );
  S *= h/3.0;

  return S;
}

double func(double x){
  return 4*sqrt(1-x*x);
  //return 8*x*x*sqrt(2-x*x);
}
