#include <stdio.h>
#include <math.h>
//cl /source-charset:shift-jis program4.2.c /EHsc
#define EPS pow(10.0,-8.0) /* epsilonの設定 */
#define NMAX 10 /* 最大反復回数 */

void newton(double x); /* Newton法    */
double f(double x);    /* f(x)の計算  */
double df(double x);   /* f'(x)の計算 */

int main(void){
	double x;
	printf("初期値x0を入力してください\n");
	scanf_s("%lf",&x);

	newton(x);
	return 0;
}

/* Newton法 */
void newton(double x){
	int n=0; double d;

	do{
		d=-f(x)/df(x);
		x=x+d;
		n++;
	}while(fabs(d)>EPS&&n<NMAX);

	if(n==NMAX){
		printf("答えは見つかりませんでした\n");
	}
	else{
		printf("答えはx=%fです．\n",x);
	}
}

/* 関数の定義 */
double f(double x){
	return x*x*x*x*x*x-14*x*x*x*x+49*x*x-36;
}

double df(double x){
	return 6*x*x*x*x*x-14*4*x*x*x+49*2*x;
}
