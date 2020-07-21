#include <stdio.h>
#include <math.h>
//cl /source-charset:shift-jis program4.2.c /EHsc
#define EPS pow(10.0,-8.0) /* epsilon�̐ݒ� */
#define NMAX 10 /* �ő唽���� */

void newton(double x); /* Newton�@    */
double f(double x);    /* f(x)�̌v�Z  */
double df(double x);   /* f'(x)�̌v�Z */

int main(void){
	double x;
	printf("�����lx0����͂��Ă�������\n");
	scanf_s("%lf",&x);

	newton(x);
	return 0;
}

/* Newton�@ */
void newton(double x){
	int n=0; double d;

	do{
		d=-f(x)/df(x);
		x=x+d;
		n++;
	}while(fabs(d)>EPS&&n<NMAX);

	if(n==NMAX){
		printf("�����͌�����܂���ł���\n");
	}
	else{
		printf("������x=%f�ł��D\n",x);
	}
}

/* �֐��̒�` */
double f(double x){
	return x-cos(x);
}

double df(double x){
	return 1.0+sin(x);
}
