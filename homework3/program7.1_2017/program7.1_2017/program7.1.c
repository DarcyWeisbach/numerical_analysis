#include <stdio.h>

/* �֐��̒�` */
double func1(double x);
double func2(double x);
/* ��`���� */
double trapezoidal( double a, double b, int n, double (*f)(double) );

int main(void)
{
  int n=100;

  printf("2.0/(x*x)��[1,2]�Őϕ����܂��B��������%d�ł�\n", n);
  printf("���ʂ�%20.15f�ł�\n",trapezoidal(1.0, 2.0, n, func1) );

  printf("4.0/(1+x*x)��[0,1]�Őϕ����܂��B��������%d�ł�\n", n);
  printf("���ʂ�%20.15f�ł�\n",trapezoidal(0.0, 1.0, n, func2) );
  
  return 0;
}

/* ��`���� */
double trapezoidal( double a, double b, int n, double (*f)(double) )
{
  double T, h; 
  int i;
 
  h = ( b - a ) /n ;  /* ���ݕ��̎w�� */

  /* ��`���� */
  T = ( (*f)(a) + (*f)(b) ) / 2.0; 
  for ( i = 1; i < n; i++) T += (*f)( a + i*h ); 
  T *= h;

  return T;
}

/* �֐��̒�` */
double func1(double x)
{ 
  return( 2.0/(x*x) );
}

double func2(double x)
{ 
  return( 4.0 / (1.0+x*x) );
}
