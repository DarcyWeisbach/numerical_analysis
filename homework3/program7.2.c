#include <stdio.h>

/* �֐��̒�` */
double func1(double x);
double func2(double x);
/* �V���v�\������ */
double simpson( double a, double b, int n, double (*f)(double) );

int main(void)
{
  int n=50;

  printf("2.0/(x*x)��[1,2]�Őϕ����܂��B��������%d�ł�\n", 2*n);
  printf("���ʂ�%20.15f�ł�\n",simpson(1.0, 2.0, n, func1) );

  printf("4.0/(1+x*x)��[0,1]�Őϕ����܂��B��������%d�ł�\n", 2*n);
  printf("���ʂ�%20.15f�ł�\n",simpson(0.0, 1.0, n, func2) );
  
  return 0;
}

/* �V���v�\������ */
double simpson( double a, double b, int n, double (*f)(double) )
{
  double S, h; 
  int i;
 
  h = ( b - a ) /( 2.0*n ) ;  /* ���ݕ��̎w�� */

  /* �V���v�\������ */
  S = ( (*f)(a) + (*f)(b) ) ;
  for ( i = 1; i < n; i++)
  {
    S += 4.0*(*f)( a + (2.0*i-1.0)*h ) + 2.0*(*f)( a + 2.0*i*h );
  }
  S += 4.0*(*f)( a + (2.0*n-1.0)*h );
  S *= h/3.0;

  return S;
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
